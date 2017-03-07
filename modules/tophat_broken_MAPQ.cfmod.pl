#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../source";
use CF::Constants;
use CF::Helpers;
use File::Copy qw(move);

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
	'cores' 	=> ['1', '6'],
	'memory' 	=> ['8G', '12G'],
	'modules' 	=> ['tophat', 'samtools'],
	'references'=> 'bowtie2',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Alignment can take ages. Be conservative..
		return CF::Helpers::minutes_to_timestamp ($num_files * 10 * 60);
	}
);

# Help text
my $helptext = "".("-"x15)."\n Tophat Module\n".("-"x15)."\n
TopHat is a fast splice junction mapper for RNA-Seq reads.
It aligns RNA-Seq reads to mammalian-sized genomes using the ultra
high-throughput short read aligner Bowtie2, and then analyzes the
mapping results to identify splice junctions between exons.\n
This module detects the encoding of input fastq files and whether
they're paired end or not. It can run with or without a supplied
GTF path but requires a bowtie2 path.

Note that this module uses the -g 1 parameter to only return unique
hits. However, a bug in tophat means that all MAPQ values are
incorrectly returned as 50. This has been reported to the tophat
authors, in the mean time, you can use the new tophat.cfmod module
to fix this behaviour.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
# Check that we have a genome defined
if(!defined($cf{'refs'}{'bowtie2'})){
	die "\n\n###CF Error: No bowtie2 path found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
} else {
	warn "\nAligning against bowtie2 path: $cf{refs}{bowtie2}\n\n";
}

# Use a GTF file if we have one
my $gtf = defined($cf{'refs'}{'gtf'}) ? "-G $cf{refs}{gtf}" : '';

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
my $version = `tophat --version`;
warn "---------- Tophat version information ----------\n";
warn $version;
warn "------- End of Tophat version information ------\n";
if($version =~ /TopHat v([\d\.]+)/){
  warn "###CFVERS tophat\t$1\n\n";
}

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});

# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

# Go through each single end files and run Tophat
my $output_dir;
if($se_files && scalar(@$se_files) > 0){
	foreach my $file (@$se_files){
		my $timestart = time;

		# Figure out the encoding if we don't already know
		if(!$encoding){
			($encoding) = CF::Helpers::fastq_encoding_type($file);
		}
		# Tophat only accepts --solexa-quals, --solexa1.3-quals and --integer-quals
		# Defaults to phred-33 so leave blank if that returned
		my $enc = "";
		if($encoding eq 'solexa'){
			$enc = '--solexa-quals';
		} elsif($encoding eq 'phred64'){
			$enc = '--solexa1.3-quals';
		} elsif($encoding eq 'integer'){
			$enc = '--integer-quals';
		}

		$output_dir = $file;
		$output_dir =~ s/.gz$//;
		$output_dir =~ s/.fq$//;
		$output_dir =~ s/.fastq$//;
		$output_dir =~ s/_[1-4]$//;
		$output_dir =~ s/_R[1-4]//;
		$output_dir .= '_tophat';

		my $output_fn = $output_dir."/accepted_hits.bam";

		my $command = "tophat -p $cf{cores} -g 1 $enc $gtf -o $output_dir $cf{refs}{bowtie2} $file";
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			# Tophat worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF Tophat (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				$output_fn = clean_output($output_dir, $file);
				print RUN "$cf{job_id}\t$output_fn\n";
			} else {
				warn "\n###CF Error! Tophat output file $output_fn not found..\n";
			}
		} else {
			warn "\n###CF Error! Tophat (SE mode) failed for input file '$file': $? $!\n\n";
		}

	}
}

# Go through the paired end files and run Tophat
if($pe_files && scalar(@$pe_files) > 0){
	foreach my $files_ref (@$pe_files){
		my @files = @$files_ref;
		if(scalar(@files) == 2){
			my $timestart = time;

			# Figure out the encoding if we don't already know
			if(!$encoding){
				($encoding) = CF::Helpers::fastq_encoding_type($files[0]);
			}
			# Tophat only accepts --solexa-quals, --solexa1.3-quals and --integer-quals
			# Defaults to phred-33 so leave blank if that returned
			my $enc = "";
			if($encoding eq 'solexa'){
				$enc = '--solexa-quals';
			} elsif($encoding eq 'phred64'){
				$enc = '--solexa1.3-quals';
			} elsif($encoding eq 'integer'){
				$enc = '--integer-quals';
			}

			$output_dir = $files[0];
			$output_dir =~ s/.gz$//;
			$output_dir =~ s/.fq$//;
			$output_dir =~ s/.fastq$//;
			$output_dir =~ s/_[1-4]$//;
			$output_dir =~ s/_R[1-4]//;
			$output_dir .= '_tophat';

			my $output_fn = $output_dir."/accepted_hits.bam";

			my $command = "tophat -p $cf{cores} -g 1 $enc $gtf -o $output_dir $cf{refs}{bowtie2} $files[0] $files[1]";
			warn "\n###CFCMD $command\n\n";

			if(!system ($command)){
				# Tophat worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF Tophat (PE mode) successfully exited, took $duration....\n";
				if(-e $output_fn){
					$output_fn = clean_output($output_dir, $files[0]);
					print RUN "$cf{job_id}\t$output_fn\n";
				} else {
					warn "\n###CF Error! Tophat output file $output_fn not found\n";
				}
			} else {
				warn "\n###CF Error! Tophat (PE mode) failed for input file '$files[0]': $? $!\n\n";
			}
		} else {
			warn "\n###CF Error! Tophat paired end files had ".scalar(@files)." input files instead of 2\n";
		}
	}
}


# Clear up Tophat output
sub clean_output {

	my ($output_dir, $file) = @_;
	my $results_fn = "$output_dir/accepted_hits.bam";

	if(-e "$output_dir/accepted_hits.bam" && -e "$output_dir/align_summary.txt"){
		if(move ("$output_dir/accepted_hits.bam", $file."_tophat.bam")){
			warn "Moved tophat output to ${file}_tophat.bam\n";
			$results_fn = $file."_tophat.bam";
			if(move ("$output_dir/align_summary.txt", $file."_align_summary.txt")){
				warn "Moved tophat summary to ${file}_align_summary.txt\n";
				if(unlink glob "$output_dir/logs/*"){
					warn "Deleted other tophat output files: $output_dir/logs/*\n";
					if(rmdir "$output_dir/logs"){
						warn "Deleted tophat logs directory $output_dir/logs\n";
						if(unlink glob "$output_dir/*"){
							warn "Deleted other tophat output files: $output_dir/*\n";
							if(rmdir $output_dir){
								warn "Deleted tophat output directory $output_dir\n";
							} else {
								warn "Could not delete tophat output directory $output_dir: $!\n";
							}
						} else {
							warn "Could not delete other tophat output files: $output_dir/*: $!\n";
						}
					} else {
						warn "Could not delete tophat logs directory $output_dir/logs: $!\n";
					}
				} else {
					warn "Could not delete other tophat log files: $output_dir/logs/*: $!\n";
				}
			} else {
				warn "Could not move tophat summary to ${file}_align_summary.txt: $!\n";
			}
		} else {
			warn "Could not move tophat output to ${file}_tophat.bam: $!\n";
		}
	} else {
		warn "Could not find tophat output files $output_dir/accepted_hits.bam and $output_dir/align_summary.txt\n";
	}

	return $results_fn;
}


close (RUN);
