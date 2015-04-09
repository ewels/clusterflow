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
	'cores' 	=> '3',
	'memory' 	=> '3G',
	'modules' 	=> ['fastqc','cutadapt','trim_galore'],
	'time' 		=> sub {
		my $runfile = $_[0];
		my $num_files = $runfile->{'num_starting_merged_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Trim Galore typically takes less than 3 hours per file
		return CF::Helpers::minutes_to_timestamp ($num_files * 4 * 60);
	}
);

# Help text
my $helptext = "".("-"x20)."\n Trim Galore Module\n".("-"x20)."\n
Trim Galore is a wrapper tool around Cutadapt and FastQC
to consistently apply quality and adapter trimming to FastQ files.\n
This module intelligently works on single end and paired end input files
and tries to calculate input file fastq encoding.\n
By default, it will not trim files with reads shorter than 50bp
(it checks the first file). This behaviour can be overwritten
with the min_readlength pipeline parameter. For example:
#trim_galore min_readlength=0\n\nThe module will use the default
Illumina common adapter sequence to trim with, but this can be
overridden by specifying adapter=[Sequence] in the options.\n\n";

# Setup
my %runfile = CF::Helpers::module_start(\@ARGV, \%requirements, $helptext);

# MODULE
my $timestart = time;

open (RUN,'>>',$runfile{'run_fn'}) or die "###CF Error: Can't write to $runfile{run_fn}: $!";

# Print version information about the module.
warn "---------- Cutadapt version information ----------\n";
warn `cutadapt --version`;
warn "\n------- End of Cutadapt version information ------\n";
warn "---------- Trim Galore! version information ----------\n";
warn `trim_galore --version`;
warn "\n------- End of Trim Galore! version information ------\n";

# Read options from the pipeline parameters
my $min_readlength = (defined($runfile{'params'}{'min_readlength'})) ? $runfile{'params'}{'min_readlength'} : 50;
my $q_cutoff = (defined($runfile{'params'}{'q'})) ? "-q ".$runfile{'params'}{'q_cutoff'} : '';
my $stringency = (defined($runfile{'params'}{'stringency'})) ? "--stringency ".$runfile{'params'}{'stringency'} : '';
my $adapter = (defined($runfile{'params'}{'adapter'})) ? "--adapter ".uc($runfile{'params'}{'adapter'}) : '';
my $RRBS = (defined($runfile{'params'}{'RRBS'})) ? "--RRBS" : '';

my $clip_r1 = "";
my $clip_r2 = "";
if(defined($runfile{'params'}{'pbat'})){
	$clip_r1 = "--clip_r1 4";
	$clip_r2 = "--clip_r2 4";
}
if(defined($runfile{'params'}{'single_cell'})){
	$clip_r1 = "--clip_r1 9";
	$clip_r2 = "--clip_r2 9";
}


# Are these reads long enough to trim? Check first file and assume same for rest
if(!CF::Helpers::fastq_min_length($runfile{'prev_job_files'}[0], $min_readlength)){
	# First file didn't have long enough reads for trimming
	# Print output filenames and exit
	foreach my $file (@{$runfile{'prev_job_files'}}){
		print RUN "$runfile{job_id}\t$file\n";
	}
	close (RUN);
	warn "###CF Trim galore didn't run as reads were too short..\n";
	exit;
}

# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%runfile, @{$runfile{'prev_job_files'}});

# Go through each single end files and run trim galore
if($se_files && scalar(@$se_files) > 0){
	foreach my $file (@$se_files){

		# Figure out the encoding if we don't already know
		if(!$encoding){
			($encoding) = CF::Helpers::fastq_encoding_type($file);
		}
		my $enc = "";
		if($encoding eq 'phred33' || $encoding eq 'phred64'){
			$enc = '--'.$encoding;
		}

		my $output_fn = trim_galore_basename($file).'_trimmed.fq.gz';

		my $fqc = (defined($runfile{'params'}{'nofastqc'})) ? '' : '--fastqc_args "-q"';

		my $command = "trim_galore --gzip $enc $stringency $RRBS $adapter $q_cutoff $clip_r1 $fqc $file";

		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			# Trim Galore worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF Trim galore (SE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				print RUN "$runfile{job_id}\t$output_fn\n";
			} else {
				warn "\n###CF Error! Trim Galore output file $output_fn not found..\n";
			}
		} else {
			warn "###CF Error! Trim Galore (SE mode) exited with an error state for input file '$file': $? $!\n\n";
		}
	}
}

# Go through the paired end files and run trim galore
if($pe_files && scalar(@$pe_files) > 0){
	foreach my $files_ref (@$pe_files){
		my @files = @$files_ref;
		if(scalar(@files) == 2){

			# Figure out the encoding if we don't already know
			if(!$encoding){
				($encoding) = CF::Helpers::fastq_encoding_type($files[0]);
			}
			if(!$encoding){
				($encoding) = CF::Helpers::fastq_encoding_type($files[1]);
			}

			my $enc = "";
			if($encoding eq 'phred33' || $encoding eq 'phred64'){
				$enc = '--'.$encoding;
			}

			my $output_fn_1 = trim_galore_basename($files[0]).'_val_1.fq.gz';
			my $output_fn_2 = trim_galore_basename($files[1]).'_val_2.fq.gz';

			my $fqc = (defined($runfile{'params'}{'nofastqc'})) ? '' : '--fastqc';

			my $command = "trim_galore --paired --gzip $enc $stringency $RRBS $adapter $q_cutoff $clip_r1 $clip_r2 $fqc $files[0] $files[1]";

			warn "\n###CFCMD $command\n\n";

			if(!system ($command)){
				# Trim Galore worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF Trim galore (PE mode) successfully exited, took $duration..\n";
				if(-e $output_fn_1){
					print RUN "$runfile{job_id}\t$output_fn_1\n";
				} else {
					warn "\n###CF Error! Trim Galore output file $output_fn_1 not found..\n";
				}
				if(-e $output_fn_2){
					print RUN "$runfile{job_id}\t$output_fn_2\n";
				} else {
					warn "\n###CF Error! Trim Galore output file $output_fn_2 not found..\n";
				}
			} else {
				warn "###CF Error! Trim Galore (PE mode) exited with an error state for input file '".$files[0]."': $? $!\n\n";
			}
		} else {
			warn "\n###CF Error! Trim Galore paired end files had ".scalar(@files)." input files instead of 2..\n";
		}
	}
}

sub trim_galore_basename {

	my ($fn) = @_;

	if ($fn =~ /\.fastq$/){
		$fn =~ s/\.fastq$//;
	} elsif ($fn =~ /\.fastq\.gz$/){
		$fn =~ s/\.fastq\.gz$//;
	} elsif ($fn =~ /\.fq$/){
		$fn =~ s/\.fq$//;
	} elsif ($fn =~ /\.fq\.gz$/){
		$fn =~ s/\.fq\.gz$//;
	}

	return $fn;
}


close (RUN);
