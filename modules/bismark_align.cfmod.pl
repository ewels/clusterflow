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
	'cores' 	=> ['4', '6'],
	'memory' 	=> ['18G', '25G'],
	'modules' 	=> ['bowtie','bowtie2','bismark','samtools'],
	'references'=> 'bismark',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Bismark alignment typically takes around than 3 hours per BAM file
		# May need tweaking. Could also drop the estimate if we're multi-threading?
		# Update: Keeps timing out. Bumping to 48 hours per file.
		return CF::Helpers::minutes_to_timestamp ($num_files * 48 * 60);
	}
);

# Help text
my $helptext = "".("-"x22)."\n Bismark Align Module\n".("-"x22)."\n
The bismark_align module runs the main bismark script.
Bismark is a program to map bisulfite treated sequencing reads
to a genome of interest and perform methylation calls.\n
PBAT, single cell and Bowtie 1 mode can be specified in
pipelines with the pbat, single_cell bt1 parameters. For example:
  #bismark_align	pbat\n
  #bismark_align	bt1\n
  #bismark_align	single_cell\n
  #bismark_align	unmapped\n
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
warn "Allocated $cf{cores} cores and $cf{memory} memory.\n";

# Read options from the pipeline parameters
my $pbat = defined($cf{'params'}{'pbat'}) ? "--pbat" : '';
my $unmapped = defined($cf{'params'}{'unmapped'}) ? "--unmapped" : '';
my $bowtie = defined($cf{'params'}{'bt1'}) ? "--bowtie1" : "--bowtie2";
my $non_directional = defined($cf{'params'}{'single_cell'}) ? "--non_directional" : '';
my $subsample = defined($cf{'params'}{'subsample'}) ? "-u 1000000" : '';

if(defined($cf{'params'}{'subsample'})){
	warn "WARNING! Bismark running in subsample mode - only first 1000000 reads will be aligned.\n";
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
		
		my $output_fn = $file;
		$output_fn =~ s/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//; # attempting to remove fastq.gz etc to make filename a little shorter 05 02 2016. Felix
		$output_fn .= "_".$cf{config}{genome};
		my $basename = $output_fn;
		if($bowtie == '--bowtie2'){
		    $output_fn .= "_bismark_bt2.bam";
		    $basename .= "_bismark_bt2";
		} else {
		    $output_fn .= "_bismark.bam";
		    $basename .= "_bismark";
		}
		
		my $command = "bismark --bam --basename $basename $bowtie $pbat $unmapped $non_directional $enc $cf{refs}{bismark} $file";
		
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
			
			my $output_fn = $files[0];
			$output_fn =~ s/(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$//; # attempting to remove fastq.gz etc to make filename a little shorter 05 02 2016. Felix
			$output_fn .= "_".$cf{config}{genome};
			my $basename = $output_fn;

			if($bowtie == '--bowtie2'){
			    $output_fn .= "_bismark_bt2_pe.bam";
			    $basename .= "_bismark_bt2_pe";
			} else {
			    $output_fn .= "_bismark_pe.bam";
			    $basename .= "_bismark_pe";
			}

			my $command = "bismark --bam --basename $basename $bowtie $pbat $unmapped $non_directional $enc $cf{refs}{bismark} -1 ".$files[0]." -2 ".$files[1];
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
