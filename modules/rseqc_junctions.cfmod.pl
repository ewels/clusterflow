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
Module to run RSeQC 'junctions annotation' and 'junction saturation' scripts:
http://rseqc.sourceforge.net/#junction-annotation-py
http://rseqc.sourceforge.net/#junction-saturation-py\n
Use the 'keep_intermediate' parameter to keep the r plotting files, otherwise
these will be deleted.\n";

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
warn `junction_annotation.py --version`;
warn `junction_saturation.py --version`;
warn "------- End of RSeQC version information ------\n";

# Set up optional parameters
my $keep_intermediate = (defined($cf{'params'}{'keep_intermediate'})) ? 1 : 0;

foreach my $file (@{$cf{'prev_job_files'}}){
	my $timestart = time;

	# Output file name prefix
	my $output_prefix = $file;
	$output_prefix =~ s/.bam//;

	# Step 1 - Junction annotation
	my $junc_annotation_cmd = "junction_annotation.py -i $file -o $output_prefix -r $cf{refs}{bed12}";
	warn "\n###CFCMD $junc_annotation_cmd\n\n";

	if(!system ($junc_annotation_cmd)){

		# Delete intermediate file
		if(!$keep_intermediate){
			unlink($output_prefix.".junction_plot.r");
		}

		my $outputfile1 = $output_prefix.".splice_events.pdf";
		my $outputfile2 = $output_prefix.".splice_junction.pdf";
		if(-e $outputfile1){
			print RUN $cf{'job_id'}."\t$outputfile1\n";
			print RUN $cf{'job_id'}."\t$outputfile2\n";
		} else {
			warn "\n###CF Error! RSeQC junction annotation output file $outputfile1 not found..\n";
		}

	} else {
		warn "\n###CF Error! RSeQC junction annotation plot failed, exited in an error state: $? $!\n\n";
	}


	# Step 2 - Junction saturation
	my $junc_saturation_cmd = "junction_saturation.py -i $file -o $output_prefix -r $cf{refs}{bed12}";
	warn "\n###CFCMD $junc_saturation_cmd\n\n";

	if(!system ($junc_saturation_cmd)){
		# command worked - print out resulting filenames
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF RSeQC junction annotation and saturation analysis finished, took $duration..\n";

		# Delete intermediate file
		if(!$keep_intermediate){
			unlink($output_prefix.".junctionSaturation_plot.r");
		}

		my $outputfile = $output_prefix.".junctionSaturation_plot.pdf";
		if(-e $outputfile){
			print RUN $cf{'job_id'}."\t$outputfile\n";
		} else {
			warn "\n###CF Error! RSeQC junction saturation plot output file $outputfile not found..\n";
		}

	} else {
		warn "\n###CF Error! RSeQC junction saturation plot failed, exited in an error state: $? $!\n\n";
	}
}


close (RUN);
