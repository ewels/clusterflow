#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../source";
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
	'memory' 	=> '10G',
	'modules' 	=> ['bowtie','samtools'],
	'references'=> 'bowtie',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Bowtie alignment typically takes less than 10 hours per BAM file
		# This is probably inaccurate? May need tweaking.
		return CF::Helpers::minutes_to_timestamp ($num_files * 14 * 60);
	}
);

# Help text
my $helptext = "".("-"x17)."\n Bowtie 1 Module\n".("-"x17)."\n
Bowtie is an ultrafast, memory-efficient short read aligner.
It aligns short DNA sequences (reads) to a reference genome. This script
works out the encoding of input files, guesses whether they're paired end
or not and runs bowtie. Output is piped through samtools to generate BAM files.\n\n".
"If input files are gzipped they will either be piped through zcat
(single end) or gunzipped, processed and zipped again (paired end).\n\n".
"To use this module to align microRNAs, use the parameter mirna, eg:
#bowtie mirna\n\n".
"For further information about bowtie, please run bowtie --help\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
# Check that we have a genome defined
if(!defined($cf{'refs'}{'bowtie'})){
    die "\n\n###CF Error: No bowtie path found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
} else {
    warn "\nAligning against ".$cf{'refs'}{'bowtie'}."\n\n";
}

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
my $version = `bowtie --version`;
warn "---------- Bowtie 1 version information ----------\n";
warn $version;
warn "\n------- End of Bowtie 1 version information ------\n";
if($version =~ /bowtie version ([\d\.]+)/){
  warn "###CFVERS bowtie\t$1\n\n";
}

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

		my $output_fn = $file;
		$output_fn =~ s/\.gz$//i;
		$output_fn =~ s/\.(fq|fastq)$//i;
		$output_fn .= "_".$cf{config}{genome};
		$output_fn .= "_bowtie.bam";

		my $gzip = "";
		if($file =~ /\.gz$/){
			# $gzip = "gzip -dc $file | ";
			$gzip = "zcat $file | ";
			$file = "-";
		}

		my $command = "$gzip bowtie -p $cf{cores} -t -m 1 $enc --strata --best -S --chunkmbs 512 $cf{refs}{bowtie} $file | samtools view -bS - > $output_fn";
		if (defined($cf{'params'}{'mirna'})) {
		    $command = "$gzip bowtie -p $cf{cores} -t -n 0 -l 15 -e 99999 -k 200 $enc --best -S --chunkmbs 512 $cf{refs}{bowtie} $file | samtools view -bS - > $output_fn"
		}
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			# Bowtie worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF Bowtie (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				print RUN "$cf{job_id}\t$output_fn\n";
			} else {
				warn "\n###CF Error! Bowtie output file $output_fn not found..\n";
			}
		} else {
			warn "\n###CF Error! Bowtie (SE mode) exited in an error state for input file '$file': $? $!\n\n";
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

			# Unzip gzipped input files
			if($files[0] =~ /.gz$/){
				warn "Unzipping input files..\n";
				system("gunzip ".$files[0]);
				system("gunzip ".$files[1]);
				$files[0] =~ s/.gz$//;
				$files[1] =~ s/.gz$//;
			}

			my $output_fn = $files[0];
			$output_fn =~ s/\.gz$//i;
			$output_fn =~ s/\.(fq|fastq)$//i;
			$output_fn .= "_".$cf{config}{genome};
			$output_fn .= "_bowtie.bam";

			my $command = "bowtie -p $cf{cores} -t -m 1 $enc --strata --best --maxins 700 -S --chunkmbs 512 $cf{refs}{bowtie} -1 ".$files[0]." -2 ".$files[1]." | samtools view -bSh 2>/dev/null - > $output_fn";

			if (defined($cf{'params'}{'mirna'})) {
			    $command = "bowtie -p $cf{cores} -t -n 0 -l 15 -e 99999 -k 200 $enc --best --maxins 700 -S --chunkmbs 512 $cf{refs}{bowtie} -1 ".$files[0]." -2 ".$files[1]." | samtools view -bSh 2>/dev/null - > $output_fn";
			}
			warn "\n###CFCMD $command\n\n";

			if(!system ($command)){
				# Bowtie worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF Bowtie (PE mode) successfully exited, took $duration..\n";
				if(-e $output_fn){
					print RUN "$cf{job_id}\t$output_fn\n";
				} else {
					warn "\n###CF Error! Bowtie output file $output_fn not found..\n";
				}
			} else {
				warn "\n###CF Error! Bowtie (PE mode) exited in an error state for input file '".$files[0]."': $? $!\n\n";
			}

			# Zip input files up again
			warn "Zipping input files..\n";
			system("gzip ".$files[0]);
			system("gzip ".$files[1]);

		} else {
			warn "\n###CF Error! Bowtie paired end files had ".scalar(@files)." input files instead of 2\n";
		}
	}
}


close (RUN);
