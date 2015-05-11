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
	'cores' 	=> ['1', '8'],
	'memory' 	=> ['3G', '80G'],
	'modules' 	=> ['preseq'],
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Preseq typically takes less than 4 hours per BAM file
		return CF::Helpers::minutes_to_timestamp ($num_files * 6 * 60);
	}
);

# Help text
my $helptext = "".("-"x17)."\n Preseq Module\n".("-"x17)."\n
Preseq is a tool to calculate the complexity of a sequencing library.
The lc_extrap mode is used. See http://smithlabresearch.org/software/preseq/
for more information.\n
Input is a sorted and indexed BAM file. Only the BAM file should be
in the run file, the .bai index file will be assumed.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# --version. Returns version information about the module.
warn "---------- Preseq version information ----------\n";
warn `preseq 2>&1 | head -n 4`;
warn "\n------- End of Preseq version information ------\n";

# Go through each file and run Preseq
foreach my $file (@{$cf{'prev_job_files'}}){
	my $timestart = time;

	# Find if PE or SE from input BAM file
    my $paired = '';
    my $mode = 'SE';
	if(CF::Helpers::is_bam_paired_end($file)){
        $paired = '-P -l 999999999999';
        $mode = 'PE';
    }

    my $output_fn = $file.".preseq";

	my $command = "preseq lc_extrap -Q -B $paired $file -o $output_fn";
	warn "\n###CFCMD $command\n\n";

	if(!system ($command)){
		# Preseq worked - print out resulting filename
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF Preseq ($mode mode) successfully exited, took $duration..\n";
		if(-e $output_fn){
			print RUN "$cf{job_id}\t$output_fn\n";
		} else {
			warn "\n###CF Error! Preseq output file $output_fn not found..\n";
		}
	} else {
		warn "\n###CF Error! Preseq ($mode mode) exited in an error state for input file '$file': $? $!\n\n";
	}
}

close (RUN);
