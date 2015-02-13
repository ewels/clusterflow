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
    print CF::Helpers::allocate_cores($required_cores, 1, 8);
    exit;
}
# --mem. Return the required memory allocation.
if($required_mem){
    print CF::Helpers::allocate_memory($required_mem, '3G', '4G');
    exit;
}
# --modules. Return csv names of any modules which should be loaded.
if($required_modules){
    print 'preseq';
    exit;
}
# --help. Print help.
if($help){
    print "".("-"x17)."\n Preseq Module\n".("-"x17)."\n
Preseq is a tool to calculate the complexity of a sequencing library.
The lc_extrap mode is used. See http://smithlabresearch.org/software/preseq/
for more information.

Input is a sorted and indexed BAM file. Only the BAM file should be
in the run file, the .bai index file will be assumed.
\n\n";
    exit;
}

# MODULE
my $timestart = time;

# Read in the input files from the run file
my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config_ref) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;

if(!defined($cores) || $cores <= 0){
    $cores = 1;
}

open (RUN,'>>',$runfile) or die "###CF Error: Can't write to $runfile: $!";

# --version. Returns version information about the module.
warn "---------- Preseq version information ----------\n";
warn `preseq 2>&1 | head -n 4`;
warn "\n------- End of Preseq version information ------\n";	

# Go through each file and run Preseq
if($files && scalar(@$files) > 0){
	foreach my $file (@$files){
		
		# Find if PE or SE from input BAM file
        my $paired = '';
		if(CF::Helpers::is_bam_paired_end($file)){
            $paired = '-P';
        }
        
        my $output_fn = $file.".preseq";
        
		my $command = "preseq lc_extrap -B $paired $file -o $output_fn";
		warn "\n###CFCMD $command\n\n";
		
		if(!system ($command)){
			# Preseq worked - print out resulting filename
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF Preseq (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				print RUN "$job_id\t$output_fn\n"; 
			} else {
				warn "\n###CF Error! Preseq output file $output_fn not found..\n";
			}
		} else {
			warn "\n###CF Error! Preseq (SE mode) exited in an error state for input file '$file': $? $!\n\n";
		}
	}
}

close (RUN);
