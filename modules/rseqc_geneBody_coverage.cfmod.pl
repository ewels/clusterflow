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
Module to run RSeQC 'Gene Body Coverage' script:
http://rseqc.sourceforge.net/#genebody-coverage-py
Calculates the average coverage across genes > 100bp long and plots
a histogram. Takes aligned BAM files as input. Requires a
BED12 reference gene model.\n
Can be run as a normal module (#rseqc_geneBody_coverage), which will
make one histogram per file. Or run as a summary module
(>rseqc_geneBody_coverage) which will plot a single histogram for
all BAM files found in the run files.\n";

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
warn `geneBody_coverage.py --version`;
warn "------- End of RSeQC version information ------\n";

# Set up optional parameters
my $keep_intermediate = (defined($cf{'params'}{'keep_intermediate'})) ? 1 : 0;

# Running as a summary module
if(defined($cf{'params'}{'summary_module'})){
	my $timestart = time;

	# Find BAMs from earlier modules
	my @bamfiles;
	while (my ($jobid, $files) = each %{$cf{'files'}}) {
		foreach my $fn (@{$files}){
			if($fn =~ /\.bam$/){
				push(@bamfiles, $fn);
			}
		}
	}

	# Output file name prefix
	my $output_prefix = $cf{'pipeline_id'};

	my $cmd .= "geneBody_coverage.py -i ".join(",", @bamfiles)." -o $output_prefix -r $cf{refs}{bed12}";
	run_command($cmd, $output_prefix);
}

# Running as a regular module
else {
	foreach my $file (@{$cf{'prev_job_files'}}){

		# Output file name prefix
		my $output_prefix = $file;
		$output_prefix =~ s/.bam//;

		my $cmd .= "geneBody_coverage.py -i $file -o $output_prefix -r $cf{refs}{bed12}";
		run_command($cmd, $output_prefix);
	}
}


close (RUN);


sub run_command {

	my $cmd = $_[0];
	my $output_prefix = $_[1];
	my $timestart = time;

	warn "\n###CFCMD $cmd\n\n";

	if(!system ($cmd)){

		# command worked - print out resulting filenames
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF RSeQC gene body coverage successfully exited, took $duration..\n";

		# Delete intermediate files
		if(!$keep_intermediate){
			# TODO - don't know if this file is generated yet...?
			# unlink($output_prefix.".geneBodyCoverage.curves.r");
		}

		my $outputfile = $output_prefix.".geneBodyCoverage.curves.pdf";
		if(-e $outputfile){
			print RUN $cf{'job_id'}."\t$outputfile\n";
		} else {
			warn "\n###CF Error! RSeQC gene body coverage output file $outputfile not found..\n";
		}

	} else {
		warn "\n###CF Error! RSeQC gene body coverage failed, exited in an error state: $? $!\n\n";
	}
}
