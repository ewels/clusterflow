#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;

##########################################################################
# Copyright 2014, Philip Ewels  (phil.ewels@babraham.ac.uk)              #
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
    'cores'     => '2',
    'memory'     => '5G',
    'modules'     => 'samtools',
    'time'         => sub {
        my $cf = $_[0];
        my $num_files = $cf->{'num_starting_merged_aligned_files'};
        $num_files = ($num_files > 0) ? $num_files : 1;
        # Deduplication typically takes less than 2 hours per BAM file
        return CF::Helpers::minutes_to_timestamp ($num_files * 120);
    }
);

# Help text
my $helptext = "".("-"x15)."\n Samtools dedup Module\n".("-"x15)."\n
This module uses the samtools package to mark and remove duplicated
sequences from mapped BAM files.  The process is actually 3 steps

1)samtools sort
2)samtools rmdup
3)samtools sort -n

The output will be deduplicated BAM files sorted by sequence name.
\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);


# MODULE
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- Samtools version information ----------\n";
warn `samtools 2>&1 | head -n 4`;
warn "------- End of Samtools version information ------\n";


# Go through each file and run Tophat
foreach my $file (@{$cf{'prev_job_files'}}){
    my $timestart = time;
            
    my $pos_sorted_fn = $file;
    $pos_sorted_fn =~ s/.bam$//i;

    my $dedup_fn = $pos_sorted_fn."_raw_dedup.bam";
    my $dedup_resorted_fn = $pos_sorted_fn."_dedup"; # Samtools automatically adds the bam in this case.
    $pos_sorted_fn .= '_sorted'; # Samtools automatically adds the bam in this case.

    my $pair_flag = CF::Helpers::is_bam_paired_end($file) ? '' : '-s';

    my $first_sort_command = "samtools sort $file $pos_sorted_fn";
    my $dedup_command = "samtools rmdup $pair_flag ${pos_sorted_fn}.bam $dedup_fn";
    my $resort_command = "samtools sort -n $dedup_fn $dedup_resorted_fn";
    
	warn "\n###CFCMD $first_sort_command\n\n";
    if(!system ($first_sort_command)){
		
        # First sort worked
        warn "\n###CFCMD $dedup_command\n\n";
        if(!system ($dedup_command)){
			
            # rmdup worked
            warn "\n###CFCMD $resort_command\n\n";
            if(!system ($resort_command)){
				
                # It all worked!
                # Clean up the intermediate files
                unlink("${pos_sorted_fn}.bam");
                unlink($dedup_fn);

                my $duration =  CF::Helpers::parse_seconds(time - $timestart);
                warn "###CF Samtools dedup successfully exited, took $duration..\n";
                print RUN "$job_id\t${dedup_resorted_fn}.bam\n"; 
            } else {
                warn "\n###CF Error! Samtools dedup Name sort failed for input file '$file': $? $!\n\n";
            }
        } else {
            warn "\n###CF Error! Samtools dedup Rmdup failed for input file '$file': $? $!\n\n";
        }
    } else {
        warn "\n###CF Error! Samtools dedup Position sort failed for input file '$file': $? $!\n\n";
    }
}


close (RUN);
