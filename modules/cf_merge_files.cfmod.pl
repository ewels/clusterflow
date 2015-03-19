#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;
use POSIX;

##########################################################################
# Copyright 2014, Philip Ewels (phil.ewels@scilifelab.se)                #
#                                                                        #
# This file is part of Cluster Flow.                                     #
#                                                                        #
# Cluster Flow is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# Cluster Flow is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with Cluster Flow.  If not, see <http://www.gnu.org/licenses/>.  #
##########################################################################

# Get Options
my $required_cores;
my $required_mem;
my $required_modules;
my $run_fn;
my $help;
my $result = GetOptions ("cores=i" => \$required_cores, "mem=s" => \$required_mem, "modules" => \$required_modules, "runfn=s" => \$run_fn, "help" => \$help);

# QSUB SETUP
# --cores i = offered cores. Return number of required cores.
if($required_cores){
	print '1';
	exit;
}
# --mem. Return the required memory allocation.
if($required_mem){
	print '3G';
	exit;
}
# --modules. Return csv names of any modules which should be loaded.
if($required_modules){
	print 'samtools';
	exit;
}
# --help. Print help.
if($help){
	print <<'HELP_TEXT';
-----------------------
 CF Merge Files Module
-----------------------
This core Cluster Flow module merges files. It can be used in two ways:

1. Setting --merge on the command line or @merge_regex in the configuration
   file. The core cf executable will use this string in a regex for each
   input file and try to group them. The module will automatically be
   prepended to your pipeline. This is useful for sequencing centres who
   split samples across lanes, for example.

2. Using the module in the middle of a pipeline, with a parameter set as
   the merge regex. For example: #cf_merge_files regex="[a-z]_(P\d+)_[12]"
   This is useful when a pipeline generates multiple files from a step
   and you would like to merge them. Note that input file lists will be
   split into parallel runs by Cluster Flow, and the module will not be
   able to merge files in parallel threads.

The module uses the regex to match filenames. A single set of parentheses
should be used for the unique grouping factor. This will be used for the
resulting output filename.

For example, consider the following regex:
[1-6]_[0-9]{6}_[a-zA-Z0-9]+_(P\d+_\d+_[12]).fastq.gz

If used with the following files:
1_120119_QOS3SBTASX_P8294_101_1.fastq.gz
2_120119_QOS3SBTASX_P8294_101_1.fastq.gz

Would result in the following merged file: P8294_101_1.fastq.gz

As such, this module can also be used to clean up file names without
actually merging anything.
HELP_TEXT

	exit;
}


# MODULE
my $timestart = time;

# Read in the input files from the run file
my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config_ref) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;

# Get regex from runfile config
my $regex;
if(defined($config{merge_regex}) && length($config{merge_regex}) > 0){
	$regex = $config{merge_regex};
}

# Try to find regex string in parameters (over-rides config)
foreach my $parameter (@$parameters){
	if(substr($parameter, 0, 5) eq 'regex'){
		$regex = substr($parameter, 6, -1);
		warn "\n\nUsing regex from pipeline parameters: $regex\n\n";
	}
}

# Check that we have a regex and that it looks ok
if(length($regex) == 0){
	die "\n\n###CF Error: No merging regex found in $runfile for job $job_id\n\n";
}
my $opening_p = () = $regex =~ /\(/g;
my $closing_p = () = $regex =~ /\)/g;
my $slashes = () = $regex =~ /[^\\]\//g;
unless($opening_p == 1 && $closing_p == 1){
	die "\n\n###CF Error: Merging regex didn't have one set of parentheses for job $job_id: $regex\n\n";
}
unless($slashes == 0){
	die "\n\n###CF Error: Merging regex shouldn't have any unescaped backslashes for job $job_id: $regex\n\n";
}

# We got this far - looking good!
warn "\n\nMerging files based on regex: $regex\n\n";



#
# Group the files by their regex match
#
my $num_starting_files = scalar(@{$files});
my %file_sets;
my @newfiles;
for my $file (@{$files}){
	my $group = 'unmatched';
	if($file =~ m/$config{merge_regex}/){
		if(defined($1)){
			$group = $1;
		}
	}
	# No match - keep file as it is and move on
	if($group eq 'unmatched'){
		push(@newfiles, $file);
		next;
	}
	# Match - add to hash
	if(!defined($file_sets{$group})){
		$file_sets{$group} = ();
	}
	push(@{$file_sets{$group}}, $file);
}


#
# Merge each set of files
#
for my $group (keys(%file_sets)) {

	# Take the new file extension from the first file
	my ($ext) = $file_sets{$group}[0] =~ /(\.[^.]+)$/;
	# get the previous extension if it's gz
	if($ext eq '.gz'){
		($ext) = $file_sets{$group}[0] =~ /(\.[^.]+\.gz)$/;
	}
	my $mergedfn = $group.$ext;
	my $command;

	# Single files (simple renaming)
	if(scalar(@{$file_sets{$group}}) == 1){
		$command = "mv ".$file_sets{$group}[0]." $mergedfn";
	}
	# BAM files
	elsif($ext =~ /\.bam/i){
		$command = "samtools merge $mergedfn ".join(' ', $file_sets{$group});
	}
	# gzip files - should be able to cat these, but some software has trouble with such files
	elsif($ext =~ /\.gz/i){
		$command = "zcat ".join(' ', $file_sets{$group})." | gzip -c > $mergedfn";
	}
	# Everything else: just cat.
	else {
		$command = "cat ".join(' ', @{$file_sets{$group}})." > $mergedfn";
	}

	warn "\n###CFCMD $command\n\n";
	if(!system ($command)){
	    # Merging worked - save the new filenames
		push(@newfiles, $mergedfn);
	} else {
	    die "\n###CF Error! Merging exited with an error state for file group $group\n, files '".join(', ', $file_sets{$group})."'\nError: $? $!\n\n";
	}
}


#
# Write new filenames to run file
#
my $merged_files = scalar(@newfiles);

open (RUN,'>>',$runfile) or die "###CF Error: Can't write to $runfile: $!";
for my $output_fn (@newfiles){
	if(-e $output_fn){
		print RUN "$job_id\t$output_fn\n";
	} else {
		warn "\n###CF Error! cf_merge_files couldn't find output file $output_fn ..\n";
	}
}
close (RUN);


#
# Finish up
#
my $duration =  CF::Helpers::parse_seconds(time - $timestart);
warn "\n###CF File merging successfully merged $num_starting_files input files into $merged_files merged files, took $duration..\n";
