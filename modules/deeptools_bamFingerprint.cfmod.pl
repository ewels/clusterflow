#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;
use File::Copy qw(move);
use POSIX;

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
        'cores'         => ['2', '4'],
        'memory'        => ['4G', '8G'],
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
my $helptext = "".("-"x15)."\n DeepTools bamFingerprint\n".("-"x15)."\n
Takes two input files: a bam file and a _crosscorrelation.txt file, 
as created by phantompeaktools. Creates a \"fingerprint\" plot of 
the read coverage, using deepTools/bamFingerprint. Outputs are are 
plot files named basename_fingerprint.png. The fragment length is set
using 1) the input parameter fragmentLength 2) the file _crosscorrelation.txt
and 3) default fragment length of 200 (if no input parameter is set 
and the file _crosscorrelation.txt doesn't exist). Currently expects
a locally installed bamFingerprint.\n\n";


# Start your engines...
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

my $defaultFragLen = 200; # If fragment length not defined, use this value. 


# Print version information about the program to be executed.
warn "---------- bamFingerprint  version information ----------\n";
warn `bamFingerprint --version`;
warn "\n------- End of bamFingerprint version information ------\n";

# Open up our run file in append mode
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Go through each file and run our command
foreach my $file (@{$cf{'prev_job_files'}}){
    # Start the clock...
    my $timestart = time;

    # Check that input is bam file
    if($file !~ /\.bam$/i){
       	warn "Skipping '$file' as not a BAM file..\n";
	next;
    }

    my $crossCorrelationFile = $file."_crosscorrelation.txt"; # name of corresponding _crosscorrelation.txt file
    
    # set fragment length:
    my $fragmentLength = -1;
    
    # If paramater fragmentLength is set, use this value
    if(defined($cf{'params'}{'fragmentLength'})){
	$fragmentLength = $cf{'params'}{'fragmentLength'};
	warn "###CF bamFingerprint: Using fragment length set as parameter: $fragmentLength.\n";
    } elsif(-e $crossCorrelationFile){ # Else, look for fragment length in _crosscorrelation.txt file 
	open(IN,"$crossCorrelationFile") or die "Cannot open $crossCorrelationFile";
	my $line1 = <IN>; # Just reads first line
	my @entries = split(/\s+/, $line1);
	my $crossCorrelationStr = $entries[2];
	my @crossCorrelations  = split(/\,/, $crossCorrelationStr);
	my $nrPeaks = @crossCorrelations;
	
	$fragmentLength = $crossCorrelations[floor($nrPeaks/2)]; ## use the middle peak
	warn "###CF bamFingerprint: Using fragment length from phantompeaktools cross correlation analysis: $fragmentLength.\n";
    } else { # Else, use default value and give warning
	$fragmentLength = $defaultFragLen;
	warn "###CF bamFingerprint: WARINING! No fragment length was set. Using default value of $fragmentLength.\n";
    }

    # Generate a nice output file name
    my $output_fn = $file."_fingerprint.png";
    
    # Run bamFingerprint, to get bigWig file
  
    my $cmd = "bamFingerprint -b $file -p $cf{cores} --fragmentLength $fragmentLength -plot $output_fn";
    warn "\n###CFCMD $cmd\n\n";
    
    
    # Try to run the command - returns 0 on success (which evaluated to false)
    if(!system ($cmd)){
	# Command worked!
	# Work out how long the processing took
	my $duration =  CF::Helpers::parse_seconds(time - $timestart);
	
	# Print a success message to the log file which will be e-mailed out
	warn "###CF bamFingerprint successfully exited, took $duration..\n";
	
	# Check we can find our output filename!
	if(-e $output_fn){
	    # Print the current job ID and the output filename to the run file
	    # This is so that subsequent modules can use this output
	    print RUN "$cf{job_id}\t$output_fn\n";
	} else {
	    # Oops - can't find the output file! Err...
	    warn "\nDeepTools bamFingerprint output file $output_fn not found..\n";
	}
    } else {
	# Command returned a non-zero result, probably went wrong...
	warn "\n###CF Error! DeepTools bamFingerprint failed for '$file': $? $!\n\n";
    }
}

# Close the run file
close (RUN);
