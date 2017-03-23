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
# Copyright 2017, Russell Hamilton (rsh46@cam.ac.uk)                     #
#   (module adapted from kallisto.cfmod.pl                               #
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
        'cores'         => '8',
        'memory'        => '15G',
        'modules'       => ['salmon', 'samtools'],
        'references'=> 'salmon',
        'time'          => sub {
                my $cf = $_[0];
                my $num_files = $cf->{'num_starting_merged_aligned_files'};
                $num_files = ($num_files > 0) ? $num_files : 1;
                return CF::Helpers::minutes_to_timestamp ($num_files * 2 * 60);
        }
);

# Help text
my $helptext = "".("-"x17)."\n Salmon\n".("-"x17)."\n
Salmon is a tool for quantifying the expression of transcripts using 
RNA-seq data. Salmon uses new algorithms (specifically, coupling the 
concept of quasi-mapping with a two-phase inference procedure) to 
provide accurate expression estimates very quickly (i.e. wicked-fast) 
and while using little memory. Salmon performs its inference using an 
expressive and realistic model of RNA-seq data that takes into 
account experimental attributes and biases commonly observed in real 
RNA-seq data.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
# Check that we have a genome defined
if(!defined($cf{'refs'}{'salmon'})){
   die "\n\n###CF Error: No salmon transcriptome index found in run file $cf{run_fn} for job $cf{job_id}. Exiting.. ###";
} else {
    warn "\nPseudoaligning against $cf{refs}{salmon}\n\n";
}

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
my $version = `salmon --version`;
warn "---------- Salmon version information ----------\n";
warn $version;
warn "\n------- End of Salmon version information ------\n";
if($version =~ /version : ([\d\.]+)/){
  warn "###CFVERS salmon\t$1\n\n";
  }

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});

# Go through each single end files and run salmon
if($se_files && scalar(@$se_files) > 0){
        foreach my $file (@$se_files){
                my $timestart = time;

                my $output_dir = $file."_salmon_output";
                my $output_fn  = $file."_salmon.bam";

                my $command = "salmon quant --libType A -p $cf{cores} -i $cf{refs}{salmon} -o $output_dir -r $file --geneMap $cf{refs}{gtf} --writeMappings | samtools view -Sb - > $output_fn";

                warn "\n###CFCMD $command\n\n";

                if(!system ($command)){
                        # salmon worked - print out resulting filenames
                        my $duration =  CF::Helpers::parse_seconds(time - $timestart);
                        warn "###CF salmon (SE mode) successfully exited, took $duration..\n";
                        if(-e $output_fn){
                                print RUN "$cf{job_id}\t$output_fn\n";
                        } else {
                                warn "\n###CF Error! salmon output file $output_fn not found..\n";
                        }
                } else {
                        warn "\n###CF Error! salmon (SE mode) failed exited in an error state for input file '$file': $? $!\n\n";
                }
        }
}

# Go through the paired end files and run salmon
if($pe_files && scalar(@$pe_files) > 0){
        foreach my $files_ref (@$pe_files){
                my @files = @$files_ref;
                if(scalar(@files) == 2){
                        my $timestart = time;

                        my $output_dir = $files[0]."_salmon_output";
                        my $output_fn  = $files[0]."_salmon.bam";

                        my $command = "salmon quant --libType A -p $cf{cores} -i $cf{refs}{salmon} -o $output_dir -1 $files[0] -2 $files[1] --geneMap $cf{refs}{gtf} --writeMappings | samtools view -Sb - > $output_fn";
                        warn "\n###CFCMD $command\n\n";

                        if(!system ($command)){
                                # salmon worked - print out resulting filenames
                                my $duration =  CF::Helpers::parse_seconds(time - $timestart);
                                warn "###CF salmon (PE mode) successfully exited, took $duration..\n";
                                if(-e $output_fn){
                                        print RUN "$cf{job_id}\t$output_fn\n";
                                } else {
                                        warn "\n###CF Error! salmon output file $output_fn not found..\n";
                                }
                        } else {
                                warn "\n###CF Error! salmon (PE mode) exited in an error state for input file '".$files[0]."': $? $!\n\n";
                        }

                } else {
                        warn "\n###CF Error! salmon paired end files had ".scalar(@files)." input files instead of 2\n";
                }
        }
}

close (RUN);
