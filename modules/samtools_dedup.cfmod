#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;

##########################################################################
# Copyright 2014, Philip Ewels  (phil.ewels@babraham.ac.uk)              #
# Copyright 2015, Simon Andrews (simon.andrews@babraham.ac.uk)           #
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
my $help;
my $result = GetOptions ("cores=i" => \$required_cores, "mem=s" => \$required_mem, "modules" => \$required_modules, "help" => \$help);

# QSUB SETUP
# --cores i = offered cores. Return number of required cores.
if($required_cores){
	print CF::Helpers::allocate_cores($required_cores, 2, 2);
	exit;
}
# --mem. Return the required memory allocation.
if($required_mem){
	print CF::Helpers::allocate_memory($required_mem, '5G', '5G');
	exit;
}
# --modules. Return csv names of any modules which should be loaded.
if($required_modules){
	print 'samtools';
	exit;
}
# --help. Print help.
if($help){
	print "".("-"x15)."\n Samtools dedup Module\n".("-"x15)."\n
This module uses the samtools package to mark and remove duplicated
sequences from mapped BAM files.  The process is actually 3 steps

1)samtools sort
2)samtools rmdup
3)samtools sort -n

The output will be deduplicated BAM files sorted by sequence name.
\n\n";
	exit;
}

# MODULE
my $timestart = time;

# Read in the input files from the run file
my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config_ref) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;


if(!defined($cores) || $cores < 1){
	$cores = 1;
}

open (RUN,'>>',$runfile) or die "Can't write to $runfile: $!";


# Go through each file and run Tophat
my $output_dir;
foreach my $file (@$files){
				
		my $pos_sorted_fn = $file;
		$pos_sorted_fn =~ s/.bam$//i;

		my $dedup_fn = $pos_sorted_fn."_raw_dedup.bam";
		my $dedup_resorted_fn = $pos_sorted_fn."_dedup";# Samtools automatically adds the bam in this case.
		$pos_sorted_fn .= '_sorted'; # Samtools automatically adds the bam in this case.

		# For the rmdup we need to specify whether the file is single or paired end to get the deduplication
		# right.  We can't rely on the single/paired flag from clusterflow since the input is a bam file, so
		# we'll test it explicitly.

		my $is_paired = test_paired_end($file);

		my $pair_flag = $is_paired ? "" : "-s";
	
		my $first_sort_command = "samtools sort $file $pos_sorted_fn";
		my $dedup_command = "samtools rmdup $pair_flag ${pos_sorted_fn}.bam $dedup_fn";
		my $resort_command = "samtools sort -n $dedup_fn $dedup_resorted_fn";
		warn "\n###CFCMD $first_sort_command\n\n";
		
		if(!system ($first_sort_command)){
		    # First sort worked
				warn "\n###CFCMD $dedup_command\n\n";
				if(!system ($dedup_command)){
						# rmdup worked
						warn "\n###CFCMD $resort_command\n\n";
						if(!system ($resort_command)){
								# It all worked

								# Clean up the intermediate files
								unlink("${pos_sorted_fn}.bam");
								unlink($dedup_fn);

								my $duration =  CF::Helpers::parse_seconds(time - $timestart);
								warn "###CF Samtools dedup successfully exited, took $duration..\n";
								print RUN "$job_id\t${dedup_resorted_fn}.bam\n"; 
						}

						else {
								warn "\n###CF Error! Name sort failed for input file '$file': $? $!\n\n";
						}

		    } 
				else {
						warn "\n###CF Error! Rmdup failed for input file '$file': $? $!\n\n";
		    }
		} 
		else {
		    warn "\n###CF Error! Position sort failed for input file '$file': $? $!\n\n";
		}
		
}


close (RUN);


sub test_paired_end  {
		my ($file) = @_;

		open (PAIRED,"samtools view $file | ") or die "Can't read $file: $!";

		my $count = 0;
		while (<PAIRED>) {
				my (undef,$flag) = split(/\t/);

				if ($flag & 0x1) {
						close PAIRED;
						return (1);
				}

				++$count;

				last if ($count >= 1000);
		}
		close PAIRED;
		return (0);


}
