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
	'memory' 	=> ['4G', '5G'],
	'modules' 	=> ['bowtie2','samtools'],
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Bowtie2 alignment typically takes less than 10 hours per BAM file
		# This is probably inaccurate? May need tweaking.
		return CF::Helpers::minutes_to_timestamp ($num_files * 14 * 60);
	}
);

# Help text
my $helptext = "".("-"x17)."\n Bowtie 2 Module\n".("-"x17)."\n
Bowtie 2 is an ultrafast and memory-efficient tool for aligning sequencing
reads to long reference sequences. It is particularly good at aligning reads
of about 50 up to 100s or 1,000s of characters, and particularly good at
aligning to relatively long (e.g. mammalian) genomes.
This script works out the encoding of input files, guesses whether they're
paired end or not and runs bowtie 2. Output is piped through samtools to
generate BAM files.\n
For further information, please run bowtie2 --help\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
# Check that we have a genome defined
if(!defined($cf{'refs'}{'bowtie2'})){
   die "\n\n###CF Error: No bowtie2 reference path found in run file $cf{run_fn} for job $cf{job_id}. Exiting.. ###";
} else {
    warn "\nAligning against $cf{refs}{bowtie2}\n\n";
}

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- Bowtie 2 version information ----------\n";
warn `bowtie2 --version`;
warn "\n------- End of Bowtie 2 version information ------\n";

# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});

# Go through each single end files and run Bowtie
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

		my $output_fn = $file."_bowtie2.bam";

		my $command = "bowtie2 -p $cf{cores} -t $enc -x $cf{refs}{bowtie2} -U $file | samtools view -bS - > $output_fn";
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			# Bowtie worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF Bowtie2 (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				print RUN "$cf{job_id}\t$output_fn\n";
			} else {
				warn "\n###CF Error! Bowtie2 output file $output_fn not found..\n";
			}
		} else {
			warn "\n###CF Error! Bowtie2 (SE mode) failed exited in an error state for input file '$file': $? $!\n\n";
		}
	}
}

# Go through the paired end files and run Bowtie
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

			my $output_fn = $files[0]."_bowtie2.bam";

			my $command = "bowtie2 -p $cf{cores} -t $enc -x $cf{refs}{bowtie2} -1 ".$files[0]." -2 ".$files[1]." | samtools view -bS - > $output_fn";
			warn "\n###CFCMD $command\n\n";

			if(!system ($command)){
				# Bowtie worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF Bowtie2 (PE mode) successfully exited, took $duration..\n";
				if(-e $output_fn){
					print RUN "$cf{job_id}\t$output_fn\n";
				} else {
					warn "\n###CF Error! Bowtie2 output file $output_fn not found..\n";
				}
			} else {
				warn "\n###CF Error! Bowtie2 (PE mode) exited in an error state for input file '".$files[0]."': $? $!\n\n";
			}

		} else {
			warn "\n###CF Error! Bowtie2 paired end files had ".scalar(@files)." input files instead of 2\n";
		}
	}
}


close (RUN);
