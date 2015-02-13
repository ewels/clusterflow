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

# Get Options
my $required_cores;
my $required_mem;
my $required_modules;
my $run_fn;
my $help;
my $result = GetOptions ("cores=i" => \$required_cores, "mem=s" => \$required_mem, "modules" => \$required_modules, "runfn" => \$run_fn, "help" => \$help);

# QSUB SETUP
# --cores i = offered cores. Return number of required cores.
if($required_cores){
	print 4;
	exit;
}
# --mem. Return the required memory allocation.
if($required_mem){
	print CF::Helpers::allocate_memory($required_mem, '3G', '10G');
	exit
}
# --modules. Return csv names of any modules which should be loaded.
if($required_modules){
	print 'bismark';
	exit;
}
# --help. Print help.
if($help){
	print "".("-"x38)."\n Bismark Methylation Extractor Module\n".("-"x38)."\n
The bismark_methXtract module runs the bismark_methylation_extractor script.
This reads in a bisulfite read alignment results file produced by the Bismark bisulfite
mapper and extracts the methylation information for individual cytosines in CpG, CHG 
and CHH context.\n".
"Use bismark_methylation_extractor --help for more information.\n\n";
	exit;
}

# MODULE
my $timestart = time;

# Read in the input files from the run file
my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config_ref) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;

open (RUN,'>>',$runfile) or die "###CF Error: Can't write to $runfile: $!";

# Print version information about the module.
warn "---------- bismark_methylation_extractor version information ----------\n";
warn `bismark_methylation_extractor --version`;
warn "\n------- End of bismark_methylation_extractor version information ------\n";	

# Set buffer size
my $buffer = "";
if($mem =~ /G$/){
	($buffer = $mem) =~ s/G//;
	$buffer = "--buffer ".($buffer - 1)."G";
}

# Read any options from the pipeline parameters
my $ignore_r1 = "";
my $ignore_r2 = "--ignore_r2 2";
my $gw_cov = "";
my $capture_regions = "";
foreach my $parameter (@$parameters){
	# Genome wide coverage reports
	if($parameter eq "gwcov"){
		if(!defined($config{references}{fasta})){
			warn "\n\n###CF Error - gwcov parameter specified but no genome fasta path found..\n\n";
		} else {
			$gw_cov = '--cytosine_report --genome_folder '.$config{references}{fasta};
		}
	}
	# Capture regions
	if(substr($parameter, 0, 16) eq "capture_regions="){
		$capture_regions = substr($parameter, 16);
	}
}
 


# Go through each file and deduplicate
if($files && scalar(@$files) > 0){
	foreach my $file (@$files){
		
		# Find if PE or SE from input BAM file
		if(CF::Helpers::is_bam_paired_end($file)){
			my $command = "bismark_methylation_extractor --multi $cores $ignore_r1 $ignore_r2 --bedGraph --counts $buffer --gzip -p --no_overlap --report $gw_cov $file";
			warn "\n###CFCMD $command\n\n";
			
			# Paired End BAM file
			if(!system ($command)){
				# Bismark worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "\n###CF Bismark methylation extractor (PE mode) successfully exited, took $duration..\n";
				# Attempt to plot visualizations
				if(length($gw_cov) > 0){
					warn "\nAttempting to plot coverage visualizations..\n";
					plot_visualizations($file);
				}
				my @output_fns = find_Xtracted_fns($file);
				if(scalar(@output_fns) > 0){
					foreach(@output_fns){
						print RUN "$job_id\t$_\n"; 
					}
				} else {
					warn "\n###CF Error! No bismark meth extrator output files found for input file '$file'..\n";
				}
			} else {
				die "\n###CF Error! Bismark MethXtractor (PE mode) exited with an error state for file '$file': $? $!\n\n";
			}
			
		} else {
			
			my $command = "bismark_methylation_extractor --multi $cores $ignore_r1 --bedGraph --counts $buffer --gzip -s --report $gw_cov $file";
			warn "\n###CFCMD $command\n\n";
			
			# Single End BAM file
			if(!system ($command)){
				# Bismark worked - print out resulting filenames
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "\n###CF Bismark methylation extractor (SE mode) successfully exited, took $duration..\n";
				my @output_fns = find_Xtracted_fns($file);
				if(scalar(@output_fns) > 0){
					foreach(@output_fns){
						print RUN "$job_id\t$_\n"; 
					}
				} else {
					warn "\n###CF Error! No bismark meth extrator output files found for input file '$file'..\n";
				}
			} else {
				die "\n###CF Error! Bismark MethXtractor (SE mode) exited with an error state for input file '$file': $? $!\n\n";
			}
			
		}
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


sub plot_visualizations {
	
	my ($file) = @_;
	
	# perl script paths
	my $coverage_script = "/home/phil/scripts/visualizations/bismark/bismark_coverage_curves.pl";
	my $windows_script = "/home/phil/scripts/visualizations/bismark/bismark_window_sizes.pl";
	
	# work out the filename that we're looking for
	$file =~ s/.[bs]am$//;
	$file .= '.CpG_report.txt';
	
	my $regions = '';
	if(length($capture_regions) > 0 and -e $capture_regions){
		$regions = "--regions $capture_regions";
	} elsif(!-e $capture_regions){
		warn "###CF Error! Can't find capture regions BED file $capture_regions - running without capture information.\n";
	}
	
	# Check that the coverage report exists
	unless(-e $file){
		warn "###CF Error! Can't find Bismark coverage report: $file\nSkipping visualisations..\n\n";
	} else {
		# Plot coverage curves
		if(-e $coverage_script){
			my $cov_cmd = "perl $coverage_script $regions --stranded $file";
			warn "###CFCMD $cov_cmd\n\n";
			if(!system ($cov_cmd)){
				warn "###CF Coverage curves plot successfully created.\n";
			} else {
				warn "###CF Error! Coverage curves plotting exited with an error code for input file '$file'.\n";
			}
		} else {
			warn "###CF Warning: Couldn't find find the coverage script: $coverage_script\nSkipping..\n";
		}
		# Plot coverage curves
		if(-e $windows_script){
			my $win_cmd = "perl $windows_script $regions $file";
			warn "###CFCMD $win_cmd\n\n";
			if(!system ($win_cmd)){
				warn "###CF Window sizes plot successfully created.\n";
			} else {
				warn "###CF Error! Window size plotting exited with an error code for input file '$file'.\n";
			}
	    } else {
   			warn "###CF Warning: Couldn't find find the windows script: $windows_script\nSkipping..\n";
   		}
	}
}


close (RUN);
