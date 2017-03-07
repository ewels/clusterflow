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

################################################################################
#
#	WARNING
#
#	Fastq Screen has 8 cores set as a config option in it's config file
#	Currently uses this, can't be set on the command line
#	If given less than 8 cores, will presumably overrun and cause problems
#
#	Command line override has been requested as a feature for Fastq Screen
#	Code is half-written ready to be put in below
#
#################################################################################

# Module requirements
my %requirements = (
	'cores' 	=> '8', # [1, 8]
	'memory' 	=> ['6G', '10G'],
	'modules' 	=> ['fastq_screen', 'bowtie2'],
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# FastQC typically takes less than 30 minutes per file
		return CF::Helpers::minutes_to_timestamp ($num_files * 60);
	}
);

# Help text
my $helptext = "".("-"x21)."\n FastQ Screen Module\n".("-"x21)."\n
FastQ Screen is a quality control tool that allows you to
take a sequence dataset and search it against a set of bowtie databases.\n
For further information, please run fastq_screen --help\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);


# MODULE
# Print version information about the module.
my $version = `fastq_screen --version`;
warn "---------- fastq_screen version information ----------\n";
warn $version;
warn "\n------- End of fastq_screen version information ------\n";
if($version =~ /fastq_screen v([\d\.]+)/){
  warn "###CFVERS fastq_screen\t$1\n\n";
}

# Parameters
my $conf = defined($cf{'params'}{'fastq_screen_config'}) ? "--conf ".$cf{'params'}{'fastq_screen_config'} : '';


# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});

# Go through each single end files and run Fastq Screen
if($se_files && scalar(@$se_files) > 0){
	foreach my $file (@$se_files){
		my $timestart = time;

		# Figure out the encoding if we don't already know
		if(!$encoding){
			($encoding) = CF::Helpers::fastq_encoding_type($file);
		}
		my $enc = "";
		if($encoding eq 'phred64'){
			$enc = '--illumina1_3';
		}

		# Uses --cores $cores when this is written into Fastq Screen
		my $command = "fastq_screen $conf --subset 100000 $enc --quiet --aligner bowtie2 $file";
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF Fastq Screen (SE mode) successfully ran, took $duration\n";
		} else {
			warn "###CF Error! Fastq Screen (SE mode) Failed for input file '$file': $? $!\n";
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
			if($encoding eq 'phred64'){
				$enc = '--illumina1_3';
			}

			# Uses --cores $cores when this is written into Fastq Screen
			my $command = "fastq_screen --subset 100000 $enc --quiet --aligner bowtie2 --paired ".$files[0]." ".$files[1];
			warn "\n###CFCMD $command\n\n";

			if(!system ($command)){
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF Fastq Screen (PE mode) successfully ran, took $duration\n";
			} else {
				warn "###CF Error! Fastq Screen (PE mode) failed for input file '".$files[0]."': $? $!\n";
			}
		} else {
			warn "\n###CF Error! Fastq Screen paired end files had ".scalar(@files)." input files instead of 2\n";
		}
	}
}
