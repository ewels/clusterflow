#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../source";
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

# Module requirements
my %requirements = (
	'cores' 	=> '1',
	'memory' 	=> '4G',
	'modules' 	=> 'sratoolkit',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Who knows? If it tries six times per file, could take ages!
		return CF::Helpers::minutes_to_timestamp ($num_files * 5 * 60);
	}
);

# Help text
my $helptext = "".("-"x23)."\n SRA FastQ Dump Module\n".("-"x23)."\n
This module uses the sra toolkit fastq-dump package
to extract .fastq files from .sra input. It gzips the fastq
files once produced.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);


# MODULE
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- FastQ Dump version information ----------\n";
warn `fastq-dump --version`;
warn "\n------- End of FastQ Dump version information ------\n";

# Go through each supplied file and run fastq-dump.
foreach my $file (@{$cf{'prev_job_files'}}){

	my $timestart = time;

	my $fn_base = substr($file, 0, -4);
	my @outputfiles = ($fn_base."_1.fastq", $fn_base."_2.fastq");

	for (my $attempt = 1; $attempt < 6; $attempt++) {

		my $command = "fastq-dump --split-files ./$file";
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "\n###CF FastQ Dump successfully exited on attempt $attempt, took $duration\n";
			# FastQ Dump worked - print out resulting filenames
			foreach my $output_fn (@outputfiles){
				if(-e $output_fn){
					# Zip input files up again
					warn "gzipping dumped fastq file $output_fn..\n";
					if(!system("gzip  $output_fn")){
						$output_fn .= '.gz';
						print RUN "$cf{job_id}\t$output_fn\n";
					} else {
						warn "Error - gzipping $output_fn exited with an error state: $? $!\nPrinting ungzipped filename to run file.\n";
						print RUN "$cf{job_id}\t$output_fn\n";
					}
				} else {
					warn "\nSRA dump file $output_fn not found.. (Probably single end?)\n\n";
				}
			}
			last;

		} else {

			# FastQ Dump failed - clean up partially dumped files
			foreach my $output_fn (@outputfiles){
				if(-e $output_fn){
					unlink $output_fn or die "Could not delete $output_fn : $!";
				}
			}
			warn "###CF Error! FastQ Dump failed on attempt $attempt for input file '$file': $? $!\n";

		}
	}
}

close (RUN);
