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
	'cores' 	=> '2',
	'memory' 	=> ['2G', '4G'],
	'modules' 	=> ['hicup'], # hicup_cluster/0.0.1.clusterdev
	'time' 		=> sub {
		my $runfile = $_[0];
		my $num_files = $runfile->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Wide ranging, but let's say 10 hours per BAM?
		return CF::Helpers::minutes_to_timestamp ($num_files * 14 * 60);
	}
);

# Help text
my $helptext = "".("-"x22)."\n HiCUP Module\n".("-"x22)."\n
HiCUP is a tool for mapping and performing quality
control on Hi-C data. The module will look for existing
digest files in the genome folder and create one if it
doesn't yet exist.

Three parameters can be set in CF for this module:
restriction enzyme (default HindIII), shortest read
length (default 150bp) and longest read length
(default 800bp). For example:
#hicup	longest=600
#hicup	shortest=180
#hicup	re1=A^GATCT,BglII

Use hicup --help for further information.\n\n";

# Setup
my %runfile = CF::Helpers::module_start(\@ARGV, \%requirements, $helptext);

# MODULE

# Make sure that we have a recent enough version of HiCUP
my $hicup_version = `hicup -v`;
unless(CF::Helpers::cf_compare_version_numbers($hicup_version, '0.5')){
	warn "\n\n###CF Error! HiCUP version $hicup_version is installed, version 0.5 is required for Cluster Flow. Exiting..\n\n";
	exit;
}

# Check that we have a genome defined
if(!defined($runfile{'refs'}{'fasta'})){
	die ("\n\n###CF Error! No genome path found in run file $runfile{run_fn} for job $runfile{job_id}. Exiting.. ###\n\n");
}
# Check that we have a bowtie path defined
if(!defined($runfile{'refs'}{'bowtie'})){
	die ("\n\n###CF Error! No bowtie path found in run file $runfile{run_fn} for job $runfile{job_id}. Exiting.. ###\n\n");
}
warn "\nUsing the following references:\n  $runfile{refs}{fasta}\n  $runfile{refs}{bowtie}\n\n";

open (RUN,'>>',$runfile{'run_fn'}) or die "###CF Error: Can't write to $runfile{run_fn}: $!";

# Print version information about the module.
warn "---------- HiCUP version information ----------\n";
warn `hicup --version`;
warn "\n------- End of HiCUP version information ------\n";

##############################
# Pipeline Parameters
##############################
my $longest = defined($runfile{'params'}{'longest'}) ? $runfile{'params'}{'longest'} : 800;
my $shortest = defined($runfile{'params'}{'shortest'}) ? $runfile{'params'}{'shortest'} : 150;
my $re1 = defined($runfile{'params'}{'re1'}) ? $runfile{'params'}{'re1'} : "A^AGCTT,HindIII";
my @re1_parts = split(",", $re1);


##############################
# Make Output Directories
#############################
# Figure out directory names
my %dirs;
my $digest_dir;

foreach my $file (@{$runfile{'prev_job_files'}}){
	my $dir = $file;
	$dir =~ s/(\.)?(r_[12])?(_val_[12])?\.f(ast)?q.+//i;
	$dirs{$dir} = 1;
}
foreach my $dir (keys (%dirs)){
	if (!-d $dir) {
		mkdir $dir;
	}
	if(!$digest_dir){
		$digest_dir = $dir;
	}
}


##############################
# Sort out the digest file
##############################
my $digest;
my @digest_files = glob($runfile{'refs'}{'fasta'}."Digest_*.txt");
foreach my $digest_file (@digest_files){
	if(index($digest_file, $re1_parts[1]) > 0){
		$digest = $digest_file;
		warn "Found a digest file in the genome path..\n";
	}
}
# No digest file found - create it
if(!$digest){
	my $timestart = time;
	warn "Couldn't find a digest file..\n";

	# Move into the first directory created for this run
	chdir($digest_dir);

	# Work out the command and run it
	my @genome_path_parts = split("/", $runfile{'refs'}{'fasta'});
	my $command = "hicup_digester -1 $re1 -g ".$genome_path_parts[-2]."_".$genome_path_parts[-1]." ".$runfile{'refs'}{'fasta'}."*.fa";
	warn "\nCreating new digest file.\n\n";
	warn "\n###CFCMD $command\n\n";
	if(!system ($command)){
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF HiCUP digest file successfully created, took $duration..\n";
	} else {
		die "###CF Error! HiCUP could not create new digest file. Exiting... Command: $command\n\n";
	}

	# Find the really obscure file name that has been created
	@digest_files = glob("Digest_*.txt");
	foreach my $digest_file (@digest_files){
		if(index($digest_file, $re1_parts[1])){
			$digest = "$digest_dir/$digest_file";
		}
	}

	# Move back out of the directory
	chdir("../");

}
# Still can't find the digest file. Die
if(!$digest){
	die "###CF Error! Could not find or create the HiCUP digest file. Exiting..\n\n";
}


##############################
# Loop through Files
##############################
# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%runfile, @{$runfile{'prev_job_files'}});
if($se_files && scalar(@$se_files) > 0){
	warn "\n###CF Error! HiCUP found ".scalar(@$se_files)." single-end files as input..\n";
}

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
				$enc = '--format '.$encoding.'-quals';
			}

			# Work out the output filename
			my $fn1_base = $files[0];
			$fn1_base =~ s/.gz$//;
			$fn1_base =~ s/.fq$//;
			$fn1_base =~ s/.fastq$//;
			my $fn2_base = $files[1];
			$fn2_base =~ s/.gz$//;
			$fn2_base =~ s/.fq$//;
			$fn2_base =~ s/.fastq$//;
			my $output_fn = "uniques_".$fn1_base."_trunc_".$fn2_base."_trunc.bam";

			# build command
			my $command = "hicup --zip --bowtie bowtie --digest $digest $enc --index $runfile{refs}{bowtie} --longest $longest --shortest $shortest --re1 $re1 --threads $runfile{cores} --filenames $files[0],$files[1]";
			warn "\n###CFCMD $command\n\n";

			if(!system ($command)){
				# HiCUP worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF HiCUP successfully exited, took $duration..\n";
				if(-e $output_fn){
					print RUN "$runfile{job_id}\t$output_fn\n";
				} else {
					warn "\n###CF Error! HiCUP output file $output_fn not found..\n";
				}
			} else {
				warn "###CF Error! HiCUP exited with an error state for input file '$files[0]': $? $!\n\n";
			}
		} else {
			warn "\n###CF Error! HiCUP paired end files had ".scalar(@files)." input files instead of 2..\n";
		}
	}
} else {
	warn "\n###CF Error! No input files found for HiCUP..\n";
}


close (RUN);
