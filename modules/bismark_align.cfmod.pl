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
	'cores' 	=> ['4', '64'],
	'memory' 	=> ['13G', '160G'],
	'modules' 	=> ['bowtie','bowtie2','bismark','samtools'],
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Bismark alignment typically takes around than 3 hours per BAM file
		# May need tweaking. Could also drop the estimate if we're multi-threading?
		return CF::Helpers::minutes_to_timestamp ($num_files * 5 * 60);
	}
);

# Help text
my $helptext = "".("-"x22)."\n Bismark Align Module\n".("-"x22)."\n
The bismark_align module runs the main bismark script.
Bismark is a program to map bisulfite treated sequencing reads
to a genome of interest and perform methylation calls.\n
PBAT, single cell and Bowtie 1/2 modes can be specified in
pipelines with the pbat, single_cell bt1 and bt2 parameters. For example:
#bismark_align	pbat2
#bismark_align	bt1\n
#bismark_align	bt2\n
#bismark_align	single_cell\n
Use bismark --help for further information.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# Extra log info for requirement refinement
foreach my $file (@{$cf{'starting_files'}}){
	my $filesize = CF::Helpers::bytes_to_human_readable(-s $file);
	warn "\n###CF Starting filesize $filesize: $file\n";
}

# MODULE

# Check that we have a genome defined
if(!defined($cf{'refs'}{'bismark'})){
   die "\n\n###CF Error: No genome bismark path found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
} else {
    warn "\nAligning against $cf{refs}{bismark}\n\n";
}

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- Bismark version information ----------\n";
warn `bismark --version`;
warn "\n------- End of Bismark version information ------\n";

# Work out how to parallelise bismark
warn "Allocated $cf{cores} cores and $cf{memory} memory.\n";
my $multi_cores = $cf{'cores'} / 8;
my $multi_mem = CF::Helpers::human_readable_to_bytes($cf{'memory'})/ CF::Helpers::human_readable_to_bytes('20G');
my $multi;
$multi = int($multi_cores) if ($multi_cores > $multi_mem); # Really Perl, no min()?
$multi = int($multi_mem) if ($multi_mem > $multi_cores);
my $multicore;
if($multi > 1){
	$multicore = "--multicore $multi";
	warn "Multi-threading alignment with $multicore\n";
} else {
	warn "Running in regular non multi-threaded mode.\n";
}


# Read options from the pipeline parameters
my $bt1 = defined($cf{'params'}{'bt1'}) ? 1 : 0;
my $bt2 = defined($cf{'params'}{'bt2'}) ? "--bowtie2" : '';
my $pbat = defined($cf{'params'}{'pbat'}) ? "--pbat" : '';
my $non_directional = defined($cf{'params'}{'single_cell'}) ? "--non_directional" : '';

# Work out whether we should use bowtie 1 or 2 by read length
if(!$bt1 && !$bt2){
	if(!CF::Helpers::fastq_min_length($cf{'prev_job_files'}[0], 75)){
		warn "First file has reads < 75bp long. Using bowtie 1 for aligning with bismark.\n";
		$bt1 = 1;
		$bt2 = "";
	} else {
		warn "First file has reads >= 75bp long. Using bowtie 2 for aligning with bismark.\n";
		$bt1 = 0;
		$bt2 = "--bowtie2";
	}
}

# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;


# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});


# Go through each single end files and run Bismark
if($se_files && scalar(@$se_files) > 0){
	foreach my $file (@$se_files){

		my $timestart = time;

		# Figure out the encoding if we don't already know
		if(!$encoding){
			($encoding) = CF::Helpers::fastq_encoding_type($file);
		}
		my $enc = "";
		if($encoding eq 'phred33' || $encoding eq 'phred64' || $encoding eq 'solexa'){
			$enc = '--'.$encoding.'-quals';
		}

		my $output_fn;
		if($bt2){
			$output_fn = $file."_bismark_bt2.bam";
		} else {
			$output_fn = $file."_bismark.bam";
		}

		my $command = "bismark $multicore --bam $bt2 $pbat $non_directional $enc $cf{refs}{bismark} $file";
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			# Bismark worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "\n###CF Bismark (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				print RUN "$cf{job_id}\t$output_fn\n";
			} else {
				warn "\n###CF Error! Bismark output file $output_fn not found..\n";
			}
		} else {
			die "\n###CF Error! Bismark alignment (SE mode) exited with an error state for file '$file': $? $!\n\n";
		}
	}
}

# Go through the paired end files and run Bismark
if($pe_files && scalar(@$pe_files) > 0){
	foreach my $files_ref (@$pe_files){
		my @files = @$files_ref;
		if(scalar(@files) == 2){

			my $timestart = time;

			# Figure out the encoding if we don't already know
			if(!$encoding){
				($encoding) = CF::Helpers::fastq_encoding_type($files[0]);
			}
			my $enc = "";
			if($encoding eq 'phred33' || $encoding eq 'phred64' || $encoding eq 'solexa'){
				$enc = '--'.$encoding.'-quals';
			}

			my $output_fn;
			if(length($bt2) > 0){
				$output_fn = $files[0]."_bismark_bt2_pe.bam";
			} else {
				$output_fn = $files[0]."_bismark_pe.bam";
			}

			my $command = "bismark $multicore --bam $bt2 $pbat $non_directional $enc $cf{refs}{bismark} -1 ".$files[0]." -2 ".$files[1];
			warn "\n###CFCMD $command\n\n";

			if(!system ($command)){
				# Bismark worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "\n###CF Bismark (PE mode) successfully exited, took $duration..\n";
				if(-e $output_fn){
					print RUN "$cf{job_id}\t$output_fn\n";
				} else {
					warn "\n###CF Error! Bismark output file $output_fn not found..\n";
				}
			} else {
				die "\n###CF Error! Bismark alignment (PE mode) exited with an error state for file '".$files[0]."': $? $!\n\n";
			}
		} else {
			warn "\n###CF Error! Bismark paired end files had ".scalar(@files)." input files instead of 2..\n";
		}
	}
}


close (RUN);
