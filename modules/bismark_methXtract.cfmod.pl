#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
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
	'cores' 	=> '4',
	'memory' 	=> ['5G', '10G'],
	'modules' 	=> ['bismark','samtools'],
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Bismark methylation extraction typically takes less than 4 hours per BAM file
		return CF::Helpers::minutes_to_timestamp ($num_files * 6 * 60);
	}
);

# Help text
my $helptext = "".("-"x38)."\n Bismark Methylation Extractor Module\n".("-"x38)."\n
The bismark_methXtract module runs the bismark_methylation_extractor script.
This reads in a bisulfite read alignment results file produced by the Bismark bisulfite
mapper and extracts the methylation information for individual cytosines in CpG, CHG
and CHH context.\n".
"Use bismark_methylation_extractor --help for more information.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- bismark_methylation_extractor version information ----------\n";
warn `bismark_methylation_extractor --version`;
warn "\n------- End of bismark_methylation_extractor version information ------\n";

# Adjust maximum buffer memory
my $alloc_mem = CF::Helpers::human_readable_to_bytes($cf{'memory'});
my $overhead_mem = CF::Helpers::human_readable_to_bytes('2G');
my $buffer_mem = CF::Helpers::bytes_to_human_readable($alloc_mem - $overhead_mem);

# Go through each file and deduplicate
foreach my $file (@{$cf{'prev_job_files'}}){
	my $timestart = time;
	
	# Is this PE? Should we try to trim r2?
	my $ignore_r2 = '';
	if(CF::Helpers::is_bam_paired_end($file)){
		$ignore_r2 = '--ignore_r2 2';
	}

	my $command = "bismark_methylation_extractor $ignore_r2 --multi $cf{cores} --buffer_size $buffer_mem --bedGraph --gzip --report $file";
	warn "\n###CFCMD $command\n\n";

	if(!system ($command)){
		# Bismark worked - print out resulting filenames
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "\n###CF Bismark methylation extractor successfully exited, took $duration..\n";
		my @output_fns = find_Xtracted_fns($file);
		if(scalar(@output_fns) > 0){
			foreach(@output_fns){
				print RUN "$cf{job_id}\t$_\n";
			}
		} else {
			warn "\n###CF Error! No bismark meth extrator output files found for input file '$file'..\n";
		}
	} else {
		die "\n###CF Error! Bismark MethXtractor exited with an error state for file '$file': $? $!\n\n";
	}
}

sub find_Xtracted_fns {

	my ($file) = @_;

	# strip BAM / SAM file extension
	$file =~ s/.[bs]am$//;

	my @files;
	# loop through all files in directory looking for matching filenames
	foreach (glob('*.txt *.gz')) {
		# \Q .. \E provides $file as a literal
		if(/^(CpG|Non_CpG|CHG|CHH)_[COBT]+_\Q$file\E\.txt(\.gz)?$/){
			push @files, $_;
		}
	}
	return @files;
}


close (RUN);
