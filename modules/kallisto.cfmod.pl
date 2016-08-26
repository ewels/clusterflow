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
# Copyright 2016, Russell Hamilton (rsh46@cam.ac.uk)
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
	'memory' 	=> '15G',
	'modules' 	=> ['kallisto'],
	'references'=> 'kallisto',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# kallisto alignment typically takes less than 10 hours per BAM file
		# This is probably inaccurate? May need tweaking.
		return CF::Helpers::minutes_to_timestamp ($num_files * 14 * 60);
	}
);

# Help text
my $helptext = "".("-"x17)."\n Kallisto\n".("-"x17)."\n
kallisto is a program for quantifying abundances of transcripts from RNA-Seq 
data, or more generally of target sequences using high-throughput sequencing 
reads. It is based on the novel idea of pseudoalignment for rapidly 
determining the compatibility of reads with targets, without the need for 
alignment. Pseudoalignment of reads preserves the key information needed for 
quantification, and kallisto is therefore not only fast, but also as accurate 
as existing quantification tools. In fact, because the pseudoalignment 
procedure is robust to errors in the reads, in many benchmarks kallisto 
significantly outperforms existing tools.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
# Check that we have a genome defined
if(!defined($cf{'refs'}{'kallisto'})){
   die "\n\n###CF Error: No kallisto transcriptome index found in run file $cf{run_fn} for job $cf{job_id}. Exiting.. ###";
} else {
    warn "\nPseudoaligning against $cf{refs}{kallisto}\n\n";
}

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- Kallisto version information ----------\n";
warn `kallisto version`;
warn "\n------- End of Kallisto version information ------\n";

# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});

# Go through each single end files and run kallisto
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

		my $output_dir = $file."_kallisto_output";
		my $output_fn = $file;
		$output_fn =~ s/\.gz$//i;
		$output_fn =~ s/\.(fq|fastq)$//i;
		$output_fn .= "_".$cf{config}{genome};
		$output_fn .= "_kallisto.bam";
                
                my $estFragmentLength = 200;
                my $est_sd            = 20;
                
		my $command = "kallisto quant -t $cf{cores} --pseudobam --single --fragment-length=$estFragmentLength  --sd=$est_sd -i $cf{refs}{kallisto} -o $output_dir -b 100 $file | samtools view -Sb - > $output_fn";
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			# kallisto worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF kallisto (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				print RUN "$cf{job_id}\t$output_fn\n";
			} else {
				warn "\n###CF Error! kallisto output file $output_fn not found..\n";
			}
		} else {
			warn "\n###CF Error! kallisto (SE mode) failed exited in an error state for input file '$file': $? $!\n\n";
		}
	}
}

# Go through the paired end files and run kallisto
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

			my $output_dir = $files[0]."_kallisto_output";
			my $output_fn = $files[0];
			$output_fn =~ s/\.gz$//i;
			$output_fn =~ s/\.(fq|fastq)$//i;
			$output_fn .= "_".$cf{config}{genome};
			$output_fn .= "_kallisto.bam";

			my $command = "kallisto quant -t $cf{cores} --pseudobam -i $cf{refs}{kallisto} -o $output_dir -b 100 $files[0] $files[1] | samtools view -Sb - > $output_fn";
			warn "\n###CFCMD $command\n\n";

			if(!system ($command)){
				# kallisto worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF kallisto (PE mode) successfully exited, took $duration..\n";
				if(-e $output_fn){
					print RUN "$cf{job_id}\t$output_fn\n";
				} else {
					warn "\n###CF Error! kallisto output file $output_fn not found..\n";
				}
			} else {
				warn "\n###CF Error! kallisto (PE mode) exited in an error state for input file '".$files[0]."': $? $!\n\n";
			}

		} else {
			warn "\n###CF Error! kallisto paired end files had ".scalar(@files)." input files instead of 2\n";
		}
	}
}


close (RUN);
