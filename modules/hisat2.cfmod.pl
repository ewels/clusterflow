#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
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

# Get Options
my $required_cores;
my $required_mem;
my $required_modules;
my $help;
my $result = GetOptions ("cores=i" => \$required_cores, "mem=s" => \$required_mem, "modules" => \$required_modules, "help" => \$help);

# QSUB SETUP
# --cores i = offered cores. Return number of required cores.
if($required_cores){
	print CF::Helpers::allocate_cores($required_cores, 7, 7);
	exit;
}
# --mem. Return the required memory allocation.
if($required_mem){
	print CF::Helpers::allocate_memory($required_mem, '8G', '8G');
	exit;
}
# --modules. Return csv names of any modules which should be loaded.
if($required_modules){
	print 'hisat2,samtools';
	exit;
}
# --help. Print help.
if($help){
	print "".("-"x15)."\n Hisat2 Module\n".("-"x15)."\n
Hisat2 is a spliced aligner for RNA-Seq data.  It is fast and memory
efficient.
This script works out the encoding of input files, guesses whether they're
paired end or not and runs hisat2. Output is piped through samtools to
generate BAM files.\n
For further information, please run hisat2 --help\n\n";
	exit;
}

# MODULE
my $timestart = time;

# Read in the input files from the run file
my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config_ref) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;

# Check that we have a genome defined

# TODO: Define a proper hisat2 path variable, even though it will probably
# end up being the same as the bowtie2 one.

if(!defined($config{bowtie2_path})){
    if (defined($config{bowtie_path})){
		warn "\n\n### CF Error: No bowtie2 path found in run file $runfile for job $job_id.\n Exiting.. ###\n";
		
		exit;
	}  else  {

	    warn "\n\n###CF Error: No bowtie2 path found in run file $runfile for job $job_id. Exiting.. ###";
	    exit;
	}
} else {
	warn "\nAligning against ".$config{bowtie2_path}."\n\n";
}

my $gtf = '';
if(defined($config{gtf_path})){

		my $splice_file = $config{gtf_path};
		$splice_file =~ s/\.gtf$/_hisat2_splices.txt/;
		warn "\nUsing GTF path: ".$splice_file."\n\n";
		$gtf = " --known-splicesite-infile ".$splice_file;
}


if(!defined($cores) || $cores < 1){
	$cores = 1;
}

open (RUN,'>>',$runfile) or die "Can't write to $runfile: $!";

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%config, @$files);

# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

# Go through each single end files and run Hisat2
if($se_files && scalar(@$se_files) > 0){
	foreach my $file (@$se_files){
		
		# Figure out the encoding if we don't already know
		if(!$encoding){
			($encoding) = CF::Helpers::fastq_encoding_type($file);
		}
		my $enc = "";
		if($encoding eq 'phred33' || $encoding eq 'phred64' || $encoding eq 'solexa'){
			$enc = '--'.$encoding.'-quals';
		}
		
		my $output_fn = $file."_hisat2.bam";
		
		# we are currently using a very high penalty score for soft-clipping (--sp 1000,1000) because Hisat2 seems to soft-clip even when it should run in --end-to-end mode
		# we are also filtering out unmapped reads (-F 4)
		# we are also filtering non-primary alignments (-F 256)
		my $command = "hisat2 --sp 1000,1000 -p $cores -t $enc -x ".$config{bowtie2_path}." $gtf -U $file | samtools view -bS -F 4 -F 256 - > $output_fn";
		# allowing softclipping 
		# my $command = "hisat2 -p $cores -t $enc -x ".$config{bowtie2_path}." $gtf -U $file | samtools view -bS -F 4 -F 256 - > $output_fn"; 		
		warn "\n###CFCMD $command\n\n";
		
		if(!system ($command)){
			# Hisat worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF Hisat2 (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				print RUN "$job_id\t$output_fn\n"; 
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
			
			# Figure out the encoding if we don't already know
			if(!$encoding){
				($encoding) = CF::Helpers::fastq_encoding_type($files[0]);
			}
			my $enc = "";
			if($encoding eq 'phred33' || $encoding eq 'phred64' || $encoding eq 'solexa'){
				$enc = '--'.$encoding.'-quals';
			}
			
			my $output_fn = $files[0]."_hisat2.bam";
			
			#  we are currently using a very high penalty score for soft-clipping (--sp 1000,1000) because Hisat2 seems to soft-clip even when it shoudl run in --end-to-end mode
			# we are also filtering out unmapped reads (-F 4), or reads where the mate was unmapped (-F 8)
			# we are also filtering non-primary alignments (-F 256)
			my $command = "hisat2 --sp 1000,1000 -p $cores -t $enc -x ".$config{bowtie_path}." --no-mixed --no-discordant $gtf -1 ".$files[0]." -2 ".$files[1]." | samtools view -bS -F 4 -F 8 -F 256 - > $output_fn";
			warn "\n###CFCMD $command\n\n";
			
			if(!system ($command)){
				# Bowtie worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF Hisat2 (PE mode) successfully exited, took $duration..\n";
				if(-e $output_fn){
					print RUN "$job_id\t$output_fn\n";
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
