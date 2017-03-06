#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../source";
use CF::Constants;
use CF::Helpers;

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
	'cores' 	=> ['2', '4'],
	'memory' 	=> ['8G', '16G'],
	'modules' 	=> 'BEDTools',
	'time' 		=> sub {

		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		# Default to 1 if none were parsed for whatever reason
		$num_files = ($num_files > 0) ? $num_files : 1;
		# To be safe, we set the time limit to 1 hour per bam file.
		return CF::Helpers::minutes_to_timestamp ($num_files * 1 * 60);
	}
);

# The help text
my $helptext = "".("-"x15)."\n Bedtools intersectNeg\n".("-"x15)."\n
Takes a bam file and a bed file with regions, given by the flag blacklistFile. 
Filters the bam file, so the only reads not present in the blacklist regions 
are included in the output. Usees bedtools intesect for this.\n\n";

# Start your engines...
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# Print version information about the program to be executed.
my $version = `bedtools --version`;
warn "---------- Bedtools version information ----------\n";
warn $version;
warn "\n------- End of bedtools version information ------\n";
if($version =~ /bedtools (v[\d\.]+)/){
  warn "###CFVERS BEDTools\t$1\n\n";
}

# Check if a file with blacklist regions was provided.
my $blacklistFile;
if(defined($cf{'refs'}{'blacklist'})){
    $blacklistFile = $cf{'refs'}{'blacklist'};
    warn "###CF bedtools intersectNeg: Using blacklist file $blacklistFile.\n";
}
if(defined($cf{'params'}{'blacklistFile'})){
    $blacklistFile = $cf{'params'}{'blacklistFile'};
    warn "###CF bedtools intersectNeg: Using blacklist file $blacklistFile.\n";
}

# Open up our run file in append mode
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Go through each file and run our command
foreach my $file (@{$cf{'prev_job_files'}}){
  # Start the clock...
  my $timestart = time;

  # Check that input is bam file
  if($file !~ /\.bam$/i){
      warn "\n###CF Error! bedtools_intersectNeg failed: bam file expected, got $file\n\n";
      next;
  }
  
  # Skip step in pipeline if no blacklist file given
  if(!$blacklistFile) {
      warn "###CF bedtools intersectNeg: Warning - no blacklist file given. Skipping this step.\n";
      print RUN "$cf{job_id}\t$file\n";
      next;
  }

  # Generate a nice output file name
  (my $output_fn = $file) =~ s/\.bam//i;
  $output_fn .= "_noblacklist.bam";

  ## Run bedtools intersect to remove reads from blacklist regions
  my $cmd = "bedtools intersect -abam $file -b $blacklistFile -v > $output_fn";
  warn "\n###CFCMD $cmd\n\n";

  # Try to run the command - returns 0 on success (which evaluated to false)
  if(!system ($cmd)){
    # Command worked!
    # Work out how long the processing took
    my $duration =  CF::Helpers::parse_seconds(time - $timestart);

    # Print a success message to the log file which will be e-mailed out
    warn "###CFBedtools intersectNeg successfully exited, took $duration..\n";

    # Check we can find our output filename!
    if(-e $output_fn){
      # Print the current job ID and the output filename to the run file
      # This is so that subsequent modules can use this output
      print RUN "$cf{job_id}\t$output_fn\n";
    } else {
      # Oops - can't find the output file! Err...
      warn "\nBedtools intersectNeg output file $output_fn not found..\n";
    }
  } else {
    # Command returned a non-zero result, probably went wrong...
    warn "\n###CF Error! Bedtools intersectNeg failed for '$file': $? $!\n\n";
  }
}

# Close the run file
close (RUN);
