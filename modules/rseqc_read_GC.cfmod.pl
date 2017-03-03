#!/usr/bin/env perl
use warnings;
use strict;
use FindBin qw($RealBin);
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
	'memory' 	=> '4G',
	'modules' 	=> 'rseqc',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Sorting + indexing typically takes less than 1 hour per BAM file
		return CF::Helpers::minutes_to_timestamp ($num_files * 60);
	}
);

# Help text
my $helptext = "".("-"x30)."\n RSeQC - Inner Distance\n".("-"x30)."\n
Module to run RSeQC 'read GC' script:
http://rseqc.sourceforge.net/#read-gc-py
Calculates a histogram showing the GC content of reads in a library.\n
Use the 'keep_intermediate' parameter to keep the GC_plot.r
file, otherwise this will be deleted.\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- RSeQC version information ----------\n";
warn `read_GC.py --version`;
warn "------- End of RSeQC version information ------\n";

# Set up optional parameters
my $keep_intermediate = (defined($cf{'params'}{'keep_intermediate'})) ? 1 : 0;

foreach my $file (@{$cf{'prev_job_files'}}){
	my $timestart = time;

	# Output file name prefix
	my $output_prefix = $file;
	$output_prefix =~ s/.bam//;

	my $cmd = "read_GC.py -i $file -o $output_prefix";
	warn "\n###CFCMD $cmd\n\n";

	if(!system ($cmd)){
		# command worked - print out resulting filenames
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF RSeQC read GC successfully exited, took $duration..\n";

		# Delete intermediate file
		if(!$keep_intermediate){
			unlink($output_prefix.".GC_plot.r");
		}

		my $outputfile = $output_prefix.".GC_plot.pdf";
		if(-e $outputfile){
			print RUN $cf{'job_id'}."\t$outputfile\n";
		} else {
			warn "\n###CF Error! RSeQC GC plot output file $outputfile not found..\n";
		}

	} else {
		warn "\n###CF Error! RSeQC GC plot failed, exited in an error state: $? $!\n\n";
	}
}


close (RUN);
