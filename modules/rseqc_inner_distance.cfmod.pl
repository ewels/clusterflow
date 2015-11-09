#!/usr/bin/env perl
use warnings;
use strict;
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
	'memory' 	=> '4G',
	'modules' 	=> 'rseqc',
	'references'=> 'bed12',
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
Module to run RSeQC 'inner distance' script:
http://rseqc.sourceforge.net/#inner-distance-py
Calculates the distance between reasd for an RNA-seq library, whilst
considering introns. Takes aligned BAM files as input. Requires a
BED12 reference gene model. Uses first 1000000 reads.\n
Use the 'keep_intermediate' parameter to keep the inner_distance.txt
and inner_distance_plot.r files, otherwise these will be deleted.\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE

# Check that we have a genome defined
if(!defined($cf{'refs'}{'bed12'})){
    die "\n\n###CF Error: No BED12 reference gene model found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
}elsif(! -e $cf{'refs'}{'bed12'}){
    die "\n\n###CF Error: BED12 reference gene model file not found: $cf{refs}{bed12}";
} else {
    warn "\nUsing BED gene model: ".$cf{'refs'}{'bed12'}."\n\n";
}

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- RSeQC version information ----------\n";
warn `inner_distance.py --version`;
warn "------- End of RSeQC version information ------\n";

# Set up optional parameters
my $keep_intermediate = (defined($cf{'params'}{'keep_intermediate'})) ? 1 : 0;

foreach my $file (@{$cf{'prev_job_files'}}){
	my $timestart = time;

	# Output file name prefix
	my $output_prefix = $file;
	$output_prefix =~ s/.bam//;

	my $cmd = "inner_distance.py -i $file -o $output_prefix -r $cf{refs}{bed12}";
	warn "\n###CFCMD $cmd\n\n";

	if(!system ($cmd)){
		# command worked - print out resulting filenames
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF RSeQC inner distance successfully exited, took $duration..\n";

		# Delete intermediate files
		if(!$keep_intermediate){
			unlink($output_prefix.".inner_distance.txt");
			unlink($output_prefix.".inner_distance_plot.r");
		}

		my $outputfile = $output_prefix.".inner_distance_plot.pdf";
		if(-e $outputfile){
			print RUN $cf{'job_id'}."\t$outputfile\n";
		} else {
			warn "\n###CF Error! RSeQC inner distance output file $outputfile not found..\n";
		}

	} else {
		warn "\n###CF Error! RSeQC inner distance failed, exited in an error state: $? $!\n\n";
	}
}


close (RUN);
