#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;
use File::Copy qw(move);

##########################################################################
# Copyright 2014, Philip Ewels (phil.ewels@babraham.ac.uk)               #
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

my %requirements = (
        'cores'         => '1',
        'memory'        => '1G',
        'modules'       => '',
        'time'          => sub {

                my $cf = $_[0];
                my $num_files = $cf->{'num_starting_merged_aligned_files'};
                # Default to 1 if none were parsed for whatever reason
                $num_files = ($num_files > 0) ? $num_files : 1;
                # Ususally takes around 5 mins per file. Set to 30 mins just to be sure.
                return CF::Helpers::minutes_to_timestamp ($num_files * 1 * 30);
        }
);


# The help text
my $helptext = "".("-"x15)."\n BedToNrf\n".("-"x15)."\n
Takes bed file and computes the NRF (non-redundant fraction).
Uses code taken from https://github.com/mel-astar/mel-ngs/blob/master/mel-chipseq/chipseq-metrics/CalBedNrf.pl
Output is a txt file with the total number of reads, the number of unique reads and the  NRF. 
This file is named basename_nrf.txt.\n\n";

# Start your engines...
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# Open up our run file in append mode
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Go through each file and run our command
foreach my $file (@{$cf{'prev_job_files'}}){
    # Start the clock...
    my $timestart = time;

    # Check that input is bed file
    if(!($file =~ /\.bed$/)){
	warn "\n###CF Error! bedToNrf failed: bed file expected, got $file\n\n";
	next;
    }
    
    # Generate a nice output file name
    my $output_fn = $file."_nrf.txt";
    
    # Go through bed file and compute NRF (non redundant fraction) of reads.
    open(IN,"<$file")||die "###CF Error! bedToNrf failed: $!\n";
    my $Tcnt=0;
    my $prev="NA";
    my $lcnt=0;
    while(<IN>){
	chomp;
	my @line=split("\t",$_);
	$lcnt++;
	my $t = join("_",@line[0..2]);
	$Tcnt++ unless($t eq $prev);
	$prev=$t;
    }
    close(IN);

    # Print results to out file: nr total reads, nr unique reads, NRF
    open(OUT,"> $output_fn")||die "###CF Error! bedToNrf failed: $!\n";
    print OUT "$Tcnt\t$lcnt\t".($Tcnt/$lcnt)."\n";
    close(OUT);
    
    
    my $duration = CF::Helpers::parse_seconds(time - $timestart);
     
    # Check we can find our output filename!
    if(-e $output_fn){
	warn "###CF bedToNrf done, took $duration..\n";
	
	# Print the current job ID and the output filename to the run file
	# This is so that subsequent modules can use this output
	print RUN "$cf{job_id}\t$output_fn\n";
    } else {
	# Command returned a non-zero result, probably went wrong...
	warn "\n###CF Error! bedToNrf failed: $? $!\n\n";
	# Oops - can't find the output file! Err...
	warn "\nbedToNrf output file $output_fn not found..\n";
    }
}

# Close the run file
close (RUN);
