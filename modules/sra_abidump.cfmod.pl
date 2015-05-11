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

# Module requirements
my %requirements = (
	'cores' 	=> '1',
	'memory' 	=> '500M',
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
my $helptext = "".("-"x23)."\n SRA SOLiD Dump Module\n".("-"x23)."\n
This module uses the sra toolkit abi-dump package
to extract  csqual and csfasta files from .sra input.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- ABI Dump version information ----------\n";
warn `abi-dump --version`;
warn "\n------- End of ABI Dump version information ------\n";

# Go through each supplied file and run abi-dump.
foreach my $file (@{$cf{'prev_job_files'}}){

	my $timestart = time;

	my $fn_base = substr($file, 0, -4);
	my @outputfiles = ($fn_base."_F3.csfasta.gz", $fn_base."_F3_QV.qual.gz");

	for (my $attempt = 1; $attempt < 6; $attempt++) {

		my $command = "abi-dump --gzip ./$file";
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "\n###CF SOLiD Dump successfully exited on attempt $attempt, took $duration\n";
			# SOLiD Dump worked - print out resulting filenames
			foreach my $output_fn (@outputfiles){
				if(-e $output_fn){
					print RUN "$cf{job_id}\t$output_fn\n";
				} else {
					warn "\n###CF Error! SRA dump files $output_fn not found..\n\n";
				}
			}
			last;

		} else {

			# SOLiD Dump failed - clean up partially dumped files
			foreach my $output_fn (@outputfiles){
				if(-e $output_fn){
					unlink $output_fn or die "Could not delete $output_fn : $!";
				}
			}
			warn "###CF Error! SOLiD Dump failed on attempt $attempt for input file '$file': $? $!\n";

		}
	}
}

close (RUN);
