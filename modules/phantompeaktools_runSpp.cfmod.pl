#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
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
	'cores' 	=> '1',
	'memory' 	=> '2G',
	'modules' 	=> 'samtools',
	'time' 		=> sub {

		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		# Default to 1 if none were parsed for whatever reason
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Picard usually takes around 30 minutes per BAM file?
		return CF::Helpers::minutes_to_timestamp ($num_files * 1 * 60);
	}
);

# The help text
my $helptext = "".("-"x15)."\nPhantompeaktools run SPP".("-"x15)."\n
Takes (deduplicated) bam files as input and runs cross correlation analysis,
using the run_spp.r script from the pantompeaktools packeage. Outputs for
each bam file are basename_crosscorrelation.pdf (a cross correlation plot)
and  basename_crosscorrelation.txt (a text file with various cross correlation
statistics). For convenience the input bam file is also passed on as output
to be used by the next script in the pipeline.\n\n";

# Start your engines...
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# Find the script we will run.
my $runSppScript = 0;
foreach my $path (split(/:/, $ENV{'PATH'})){
    if(-e "$path/run_spp_nodups.R"){
        $runSppScript = "$path/run_spp_nodups.R";
        last;
    }
}
die "###CF Error - could not find phantompeaktools run_spp_nodups.R script" unless($runSppScript);

# Open up our run file in append mode
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Go through each file and run our command
foreach my $file (@{$cf{'prev_job_files'}}){
  # Start the clock...
  my $timestart = time;

  # Check that input is bam file
  if($file !~ /\.bam$/i){
      warn "\n###CF Error! Phantompeaktools failed: bam file expected, got $file\n\n";
      next;
  }

  # Generate a nice output file name
  my $outputPlot = $file."_crosscorrelation.pdf";
  my $outputStats = $file."_crosscorrelation.txt";


  # Convert to bed format
  my $cmd = "Rscript $runSppScript -c=$file -savp=$outputPlot -out=$outputStats";
  warn "\n###CFCMD $cmd\n\n";

  # Try to run the command - returns 0 on success (which evaluated to false)
  if(!system ($cmd)){
      # Command worked!
      # Work out how long the processing took
      my $duration =  CF::Helpers::parse_seconds(time - $timestart);

      # Print a success message to the log file which will be e-mailed out
      warn "###CF Phantompeaktools run SPP successfully exited, took $duration..\n";

      # Check if we can find the output plot!
      if(!(-e $outputPlot)){
	  # Oops - can't find the output file! Err...
	  warn "\n Phantompeaktools run SPP output file $outputPlot not found..\n";
      }
      # Check we can find the output stats!
      if(-e $outputStats){
	  # Print the current job ID and the output filename to the run file
	  # This is so that subsequent modules can use this output
	  print RUN "$cf{job_id}\t$outputStats\n"; ## text file with statistics
	  print RUN "$cf{job_id}\t$file\n"; ## also pass on input bam file as output

      } else {
	  # Oops - can't find the output file! Err...
	  warn "\n Phantompeaktools run SPP output file $outputStats not found..\n";
      }
  } else {
      # Command returned a non-zero result, probably went wrong...
      warn "\n###CF Error! Phantompeaktools run SPP failed for '$file': $? $!\n\n";
  }
}

# Close the run file
close (RUN);
