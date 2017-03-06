#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../source";
use CF::Constants;
use CF::Helpers;

##########################################################################
# Copyright 2014, Philip Ewels (phil.ewels@babraham.ac.uk)               #
# Copyright 2015, Simon Andrews (simon.andrews@babraham.ac.uk)           #
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
	'cores' 	=> '7',
	'memory' 	=> '8G',
	'modules' 	=> ['hisat2', 'samtools'],
	'references'=> 'hisat2',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Quicker than Tophat!
		return CF::Helpers::minutes_to_timestamp ($num_files * 4 * 60);
	}
);

# Help text
my $helptext = "".("-"x15)."\n Hisat2 Module\n".("-"x15)."\n
Hisat2 is a spliced aligner for RNA-Seq data.  It is fast and memory
efficient.
This script works out the encoding of input files, guesses whether they're
paired end or not and runs hisat2. Output is piped through samtools to
generate BAM files.\n
For further information, please run hisat2 --help\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);


# MODULE
# Check that we have a genome defined
if(!defined($cf{'refs'}{'hisat2'})){
	die "\n\n###CF Error: No hisat2 ref path found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
} else {
	warn "\nAligning against hisat2 path: $cf{refs}{hisat2}\n\n";
}

# Use a splices file if we have one
my $splices = '';
if(defined($cf{'refs'}{'hisat2_splices'}) && -e $cf{'refs'}{'hisat2_splices'}){
	$splices = " --known-splicesite-infile ".$cf{'refs'}{'hisat2_splices'};
	warn "\nUsing GTF path: ".$cf{'refs'}{'hisat2_splices'}."\n\n";
} elsif(defined($cf{'refs'}{'hisat2_splices'})){
	warn "\nWarning! Splice file reference ".$cf{'refs'}{'hisat2_splices'}." specified, but file doesn't exist. Ignoring.\n\n";
}

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- Hisat2 version information ----------\n";
warn `hisat2 --version`;
warn "\n------- End of Hisat2 version information ------\n";

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});

# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

# Go through each single end files and run Hisat2
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
		$output_fn .= "_hisat2.bam";
		
		# we are currently using a very high penalty score for soft-clipping (--sp 1000,1000) because Hisat2 seems to soft-clip even when it should run in --end-to-end mode
		# we are also filtering out unmapped reads (-F 4)
		# we are also filtering non-primary alignments (-F 256)
		my $command = "hisat2 --sp 1000,1000 -p $cf{cores} -t $enc -x $cf{refs}{hisat2} $splices -U $file | samtools view -bS -F 4 -F 256 - > $output_fn";
		# allowing softclipping 
		# my $command = "hisat2 -p $cf{cores} -t $enc -x $cf{refs}{hisat2} $splices -U $file | samtools view -bS -F 4 -F 256 - > $output_fn"; 		
		warn "\n###CFCMD $command\n\n";
		
		if(!system ($command)){
			# Hisat2 worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF Hisat2 (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				print RUN "$cf{job_id}\t$output_fn\n"; 
			} else {
				warn "\n###CF Error! Hisat2 output file $output_fn not found..\n";
			}
		} else {
			warn "\n###CF Error! Hisat2 (SE mode) failed exited in an error state for input file '$file': $? $!\n\n";
		}
	}
}

# Go through the paired end files and run Hisat2
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
			$output_fn =~ s/\.gz$//i;
			$output_fn =~ s/\.(fq|fastq)$//i;
			$output_fn .= "_".$cf{config}{genome};
			$output_fn .= "_hisat2.bam";
			
			#  we are currently using a very high penalty score for soft-clipping (--sp 1000,1000) because Hisat2 seems to soft-clip even when it shoudl run in --end-to-end mode
			# we are also filtering out unmapped reads (-F 4), or reads where the mate was unmapped (-F 8)
			# we are also filtering non-primary alignments (-F 256)
			my $command = "hisat2 --sp 1000,1000 -p $cf{cores} -t $enc -x $cf{refs}{hisat2} --no-mixed --no-discordant $splices -1 ".$files[0]." -2 ".$files[1]." | samtools view -bS -F 4 -F 8 -F 256 - > $output_fn";
			warn "\n###CFCMD $command\n\n";
			
			if(!system ($command)){
				# Hisat2 worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF Hisat2 (PE mode) successfully exited, took $duration..\n";
				if(-e $output_fn){
					print RUN "$cf{job_id}\t$output_fn\n";
				} else {
					warn "\n###CF Error! Hisat2 output file $output_fn not found..\n";
				}
			} else {
				warn "\n###CF Error! Hisat2 (PE mode) exited in an error state for input file '".$files[0]."': $? $!\n\n";
			}
			
		} else {
			warn "\n###CF Error! Hisat2 paired end files had ".scalar(@files)." input files instead of 2\n";
		}
	}
}


close (RUN);
