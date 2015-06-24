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
	'cores' 	=> '2',
	'memory' 	=> '8G',
	'modules' 	=> 'picard',
	'time' 		=> sub {

		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		# Default to 1 if none were parsed for whatever reason
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Picard usually takes around 30 minutes per bam file.
		return CF::Helpers::minutes_to_timestamp ($num_files * 1 * 60);
	}
);

# The help text
my $helptext = "".("-"x15)."\n Picard Deduplicate\n".("-"x15)."\n
Uses MarkDuplicates.jar from Picard to remove duplicates from a bam file.
Input is a sorted bam file. Output is basename_dedup.bam.
A file with metrics basename_picardDupMetrics.txt is also produced..\n\n";

# Start your engines...
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# Hack to get the correct command, which changes for different picard versions
my $p = $ENV{'PICARD_HOME'};
my $md = (-e "$p/picard.jar") ? "$p/picard.jar MarkDuplicates" : "$p/MarkDuplicates.jar";

# Print version information about the program to be executed.
warn "---------- Picard version information ----------\n";
warn `java -jar $md --version`;
warn "\n------- End of picard version information ------\n";

# Open up our run file in append mode
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Go through each file and run our command
foreach my $file (@{$cf{'prev_job_files'}}){
  # Start the clock...
  my $timestart = time;

  # Check that input is bam file
  if($file !~ /\.bam$/i){
      warn "\n###CF Error! Picard Dedup failed: bam file expected, got $file\n\n";
      next;
  }

  # Generate a nice output file name
  my $output_fn = $file."_dedup.bam";
  my $metricsFile = $file."_picardDupMetrics.txt";
  
  # Remove duplicates
  my $cmd = "java -Xmx2g -jar $md INPUT=$file OUTPUT=$output_fn ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=$metricsFile VALIDATION_STRINGENCY=LENIENT";
  warn "\n###CFCMD $cmd\n\n";
  
  # Try to run the command - returns 0 on success (which evaluated to false)
  if(!system ($cmd)){
    # Command worked!
    # Work out how long the processing took
    my $duration =  CF::Helpers::parse_seconds(time - $timestart);
    
    # Print a success message to the log file which will be e-mailed out
    warn "###CF Picard Dedup successfully exited, took $duration..\n";
    
    # Check we can find our output filename!
    if(-e $output_fn){
      # Print the current job ID and the output filename to the run file
      # This is so that subsequent modules can use this output
      print RUN "$cf{job_id}\t$output_fn\n";
    } else {
      # Oops - can't find the output file! Err...
      warn "\nPicard Dedup output file $output_fn not found..\n";
    }
  } else {
    # Command returned a non-zero result, probably went wrong...
    warn "\n###CF Error! Picard Dedup failed: $? $!\n\n";
  }
}


# Close the run file
close (RUN);
