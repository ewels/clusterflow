

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
	exit;
}

if($help){
	die "\nThis is a core module which merges FastQ files for Cluster Flow.\n";
}


# MODULE
my $timestart = time;

# Read in the input files from the run file
my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config_ref) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;

# Check that we have a merge regex defined
if(!defined($config{merge_regex})){
   die "\n\n###CF Error: No merging regex found in $runfile for job $job_id. Exiting.. ###";
} else {
	warn "\n\nMerging files based on regex: ".$config{merge_regex}."\n\n";
}

# Group the files by their regex match
my $starting_files = scalar(@{$files});
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

# Merge each set of files
for my $group (keys(%file_sets)) {
	# Take the new file extension from the first file
	# We can cat gzipped files together
	my ($ext) = $file_sets{$group}[0] =~ /(\.[^.]+)$/;
	my $mergedfn = $group.$ext;
	my $command = "cat ".join(' ', $file_sets{$group})." > $mergedfn";
	warn "\n###CFCMD $command\n\n";
	if(!system ($command)){
	    # Merging worked - save the new filenames
		push(@newfiles, $mergedfn);
	} else {
	    die "\n###CF Error! FastQ merging exited with an error state for file group $group\n, files '".join(', ', $file_sets{$group})."'\nError: $? $!\n\n";
	}
}

# Write new filenames to run file
my $merged_files = scalar(@newfiles);

open (RUN,'>>',$runfile) or die "###CF Error: Can't write to $runfile: $!";
for my $output_fn (@newfiles){
	if(-e $output_fn){
		print RUN "$job_id\t$output_fn\n";
	} else {
		warn "\n###CF Error! FastQ merge couldn't find output file $output_fn ..\n";
	}
}
close (RUN);

# Finish up
my $duration =  CF::Helpers::parse_seconds(time - $timestart);
warn "\n###CF FastQ merging successfully merged $starting_files input files into $merged_files merged files, took $duration..\n";
