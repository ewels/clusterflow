#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;

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
my $result = GetOptions ("cores=i" => \$required_cores, "mem=s" => \$required_mem, "modules" => \$required_modules, "runfn" => \$run_fn, "help" => \$help);

# QSUB SETUP
# --cores i = offered cores. Return number of required cores.
if($required_cores){
	print 1;
	exit;
}
# --mem. Return the required memory allocation.
if($required_mem){
	print '1G';
	exit;
}
# --modules. Return csv names of any modules which should be loaded.
if($required_modules){
	print 'samtools,htseq';
	exit;
}
# --help. Print help.
if($help){
	print "".("-"x15)."\n HTSeq Module\n".("-"x15)."\n
Counts reads overlapping exons in a GTF file. Takes an aligned BAM
file as input and returns an annotated BAM file plus a counts file.

Use parameter stranded or stranded_rev to count in stranded or reverse
stranded modes. Defaults to not stranded.\n\n";
	exit;
}

# MODULE
my $timestart = time;

# Read in the input files from the run file
my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config_ref) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;

# Check that we have a GTF file defined
if(!defined($config{references}{gtf})){
   die "\n\n###CF Error: No GTF path found in run file $runfile for job $job_id. Exiting.. ###";
} else {
    warn "\nUsing GTF file ".$config{references}{gtf}."\n\n";
}

open (RUN,'>>',$runfile) or die "###CF Error: Can't write to $runfile: $!";

# Print version information about the module.
warn "---------- HTSeq version information ----------\n";
warn `htseq-count --help 2>&1 | tail -n 4`;
warn "\n------- End of HTSeq version information ------\n";	

# Read any options from the pipeline parameters
my $stranded = "-s no";
foreach my $parameter (@$parameters){
	if($parameter eq "stranded"){
		$stranded = "-s yes";
	}
	if($parameter eq "stranded_rev"){
		$stranded = "-s reverse";
	}
}


# Go through each supplied file and run FastQC.
foreach my $file (@$files){
	my $annotated_file = $file."_annotated.sam";
	my $counts_file = $file."_counts.txt";
	my $command = "samtools view -h $file | htseq-count -o $annotated_file -t exon $stranded -q -i 'ID' - ".$config{references}{gtf}." | sort -n -k 2 -r > $counts_file";
	warn "\n###CFCMD $command\n\n";
	
	if(!system ($command)){
		print RUN "$job_id\t$annotated_file\n";
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF HTSeq successfully exited, took $duration\n";
		if(-e $annotated_file){
			print RUN "$job_id\t$annotated_file\n"; 
		} else {
			warn "\n###CF Error! Annotated BAM output file $annotated_file not found..\n";
		}
		if(-e $counts_file){
			print RUN "$job_id\t$counts_file\n"; 
		} else {
			warn "\n###CF Error! HTSeq counts file $counts_file not found..\n";
		}
	} else {
		print "###CF Error! Error - HTSeq Failed for input file '$file': $? $!\n";
	}
}

close (RUN);