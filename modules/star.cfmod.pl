#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;

##########################################################################
# Copyright 2014, Stuart Archer                                          #
# Derivative of: bowtie2.cfmod (Copyright 2014, Philip Ewels)            #
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
	'cores' 	=> ['4', '8'],
	'memory' 	=> sub {
		my $cf = $_[0];
		my $minmem = '32G';
		my $maxmem = '1.8T';
		if (defined($cf->{'refs'}{'star'}) && -e $cf->{'refs'}{'star'}."/SA") {
			# Make a guess about memory requirement from genome
			my $estmin = int(1.2 * -s $cf->{'refs'}{'star'}."/SA");
			$minmem = CF::Helpers::allocate_memory($estmin, '8G', $maxmem);
    	}
		return CF::Helpers::bytes_to_human_readable(CF::Helpers::allocate_memory($cf->{'memory'}, $minmem, $maxmem));
	},
	'modules' 	=> ['STAR','samtools'],
	'references'=> 'star',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Alignment typically takes less than 4 hours per BAM file
		return CF::Helpers::minutes_to_timestamp ($num_files * 5 * 60);
	}
);

# Help text
my $helptext = "".("-"x17)."\n STAR Module\n".("-"x17)."\n
STAR aligner. Requires a minimum of ~30GB RAM for human genome.\n
For further information, please run STAR --help\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);



# MODULE

# Check that we have a genome defined
if(!defined($cf{'refs'}{'star'})){
	die "\n\n###CF Error: No star path found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
} else {
	warn "\nAligning against $cf{refs}{star}\n\n";
}

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- STAR version information ----------\n";
warn `STAR --version`;
warn "\n------- End of STAR version information ------\n";

# Load parameters
my $genomeLoad = (defined($cf{'params'}{'LoadAndRemove'})) ? "LoadAndRemove" : 'NoSharedMemory';
$genomeLoad = (defined($cf{'params'}{'LoadAndKeep'})) ? "LoadAndKeep" : $genomeLoad;
my $sam_attributes = (defined($cf{'params'}{'outSAMattributes'})) ? $cf{'params'}{'outSAMattributes'} : 'Standard';

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});

# FastQ encoding type. Once found on one file will assume all others are the same
my $encoding = 0;

# Go through each single end file and run STAR
foreach my $file (@$se_files){

	my $timestart = time;

	# Figure out the encoding if we don't already know
	if(!$encoding){
		($encoding) = CF::Helpers::fastq_encoding_type($file);
	}
	my $enc = "";
	my %convert_enc = ('phred33' => '0', 'phred64' => '-31', 'solexa' => '-31');  # I *think* this is correct
	if($encoding eq 'phred33' or $encoding eq 'phred64' or $encoding eq 'solexa'){
		$enc = '--outQSconversionAdd '.$convert_enc{$encoding};
	}

	# Do we need to zcat this?
	my $zcat = '';
	if ($file =~ /\.gz$/) {
		$zcat = " --readFilesCommand zcat";
	}

	my $prefix = $file;
	$prefix =~ s/\.gz$//;
	$prefix =~ s/\.fastq$//;
	$prefix =~ s/\.fq$//;
	$prefix =~ s/\_R1_001$//;
	$prefix =~ s/\_R1$//i;
	$prefix =~ s/\_val_1$//;
	my $output_fn = $prefix."_star_aligned.bam";

	my $command = "STAR --runThreadN $cf{cores} $enc --outSAMattributes $sam_attributes --genomeLoad $genomeLoad $zcat --genomeDir $cf{refs}{star} --readFilesIn $file --outFileNamePrefix $prefix --outStd SAM | samtools view -bS - > $output_fn";
	warn "\n###CFCMD $command\n\n";

	if(!system ($command)){
		# STAR worked - print out resulting filenames
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF STAR (SE mode) successfully exited, took $duration..\n";
		if(-e $output_fn){
			print RUN "$cf{job_id}\t$output_fn\n";
		} else {
			warn "\n###CF Error! star output file $output_fn not found..\n";
		}
	} else {
		warn "\n###CF Error! star (SE mode) failed exited in an error state for input file '$file': $? $!\n\n";
	}
}

# Go through the paired end files and run STAR
foreach my $files_ref (@$pe_files){
	my @files = @$files_ref;
	if(scalar(@files) == 2){

		my $timestart = time;

		# Figure out the encoding if we don't already know
		if(!$encoding){
			($encoding) = CF::Helpers::fastq_encoding_type($files[0]);
		}
		my $enc = "";
		my %convert_enc = ('phred33' => '0', 'phred64' => '-31', 'solexa' => '-31');  # I think this is correct
		if($encoding eq 'phred33' || $encoding eq 'phred64' || $encoding eq 'solexa'){
			$enc = '--outQSconversionAdd '.$convert_enc{$encoding};
		}

		my $prefix = $files[0];
		$prefix =~ s/\.gz$//;
		$prefix =~ s/\.fastq$//;
		$prefix =~ s/\.fq$//;
		$prefix =~ s/\_R1_001$//;
		$prefix =~ s/\_R1$//i;
		$prefix =~ s/\_val_1$//;
		my $output_fn = $prefix."_star_aligned.bam";

		# Do we need to zcat this?
		my $zcat = '';
		if ($files[0] =~ /\.gz$/) {
			$zcat = " --readFilesCommand zcat";
		}

		my $command = "STAR --runThreadN $cf{cores} $enc --outSAMattributes $sam_attributes --genomeLoad $genomeLoad $zcat --genomeDir $cf{refs}{star} --readFilesIn $files[0] $files[1] --outFileNamePrefix $prefix --outStd SAM | samtools view -bS - > $output_fn";
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			# STAR worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF STAR (PE mode) successfully exited, took $duration..\n";
			if(-e $output_fn){
				print RUN "$cf{job_id}\t$output_fn\n";
			} else {
				warn "\n###CF Error! STAR output file $output_fn not found..\n";
			}
		} else {
			warn "\n###CF Error! STAR (PE mode) exited in an error state for input file '$files[0]': $? $!\n\n";
		}

	} else {
		warn "\n###CF Error! STAR paired end files had ".scalar(@files)." input files instead of 2\n";
	}
}


close (RUN);
