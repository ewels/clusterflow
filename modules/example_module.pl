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

#
# NOTE - this is an example module for Cluster Flow.
#        It does not have any real functionality and should not be used in
#        production. To avoid its accidental use in Cluster Flow, the
#        file name does not end in .cfmod
#        If developing a new module from this example, remember to change
#        the filename to MODULENAME.cfmod(.ext)
#

# Cluster Flow modules are called twice - first by the core cf program
# to decide what cluster resources the module requires, then again when
# it actually runs.

# Example requirement request:
# $ module_name.cfmod.pl --requirements --run_fn thisrunfile.run --cores 6 --mem 128G
# Note that --cores and --mem give the recommended maximum resources
# Module should return it's requested resources accordingly, eg: (without hash comments)
# cores: 4
# memory: 10G
# modules: bowtie, samtools
# time: 14:00

# Example running command:
# $ module_name.cfmod.pl --run_fn thisrunfile.run --job_id this_job_id --prev_job_id prev_job_id --cores 2 --mem 8.0G
# Here, --cores and --mem give the actual allocated resources

# The Helpers CF module handles all of this through the module_start() function,
# so you don't need to worry about any of it. Just supply the function with the
# command line arguments array (\@ARGV), a hash of requirements and help text.



# First, we set the module requirements. These should be supplied as a hash with
# four keys: cores, memory, modules and time. Each can be specified as a string
# (will be printed verbatim), an array (min and max for cores and memory) or a
# subroutine. Subroutines will supplied with a single variable - the %cf hash
# reference with all available info at the time.
my %requirements = (
	'cores' 	=> ('1', '8'), # number in this range will be chosen according to --cores
	'memory' 	=> ('1G', '10G'), # number in this range will be chosen according to --mem
	'modules' 	=> 'bowtie', # can be supplied as an array of multiple modules
	'time' 		=> sub {

		# Get the runfile hash reference
		my $cf = $_[0];

		# Example to print the available keys in the supplied hash. Typical
		# useful keys: {'num_starting_files'}, {'num_starting_merged_files'}, {'num_starting_merged_aligned_files'}
		# {'pipeline_name'}, {'config'}{'force_paired_end'}, {'refs'}{'reftype'}, {'config'}{'email'}
		use Data::Dumper;
		warn Dumper($cf);

		# This key estimates how many files there will be after merging
		# and alignment (eg. no separate paired end files)
		my $num_files = $cf->{'num_starting_merged_aligned_files'};

		# Default to 1 if none were parsed for whatever reason
		$num_files = ($num_files > 0) ? $num_files : 1;

		# Make a conservative guess about how long the module will take
		# per file. Calculate minutes and pass to the minutes_to_timestamp()
		# helper function to give a properly formatted time stamp.
		# Here we request 2 hours per BAM file for example.
		return CF::Helpers::minutes_to_timestamp ($num_files * 2 * 60);
	}
);

# The help text - this is returned if the module is called with --help
# Cluster Flow does this if cf --help module_name is called
my $helptext = "".("-"x15)."\n Example Perl Module\n".("-"x15)."\n
This is a dummy module which provides example code to be used
by anyone interested in writing their own Cluster Flow module.
It doesn't do anything, but its source code is nicely commented..\n\n";

# Pass everything to the helper function. This parses all of the information
# possible from the supplied run file. If --requirements or --help has been
# passed, it gives the relevant output and exits. If the module is running
# "for real", it returns a hash (not a hash reference) into %cf with
# lots of useful information needed for the module.
my %cf = CF::Helpers::module_start(\%requirements, $helptext);




#
# MODULE EXECUTION CODE
#

# Check that we have a reference genome defined, die if we don't.
if(!defined($cf{'refs'}{'fasta'})){
	# Any STDERR / STDOUT prefixed with ###CF will be highlighted in the summary e-mails
	die "\n\n###CF Error: No genome fasta path found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
} else {
	# This goes into the log file, but not the summary e-mail
    warn "\nUsing the genome path: ".$cf{'refs'}{'fasta'}."\n\n";
}

# Print version information about the program to be executed.
warn "---------- < module > version information ----------\n";
warn `MY_COMMAND --version`;
warn "\n------- End of < module > version information ------\n";

# Read any options from the pipeline parameters - these can come from
# pipeline files (eg #mymodule	myparameter) or from the command line
# at run time (eg. cf --param myparamter mymodule *gz)
# Typically, this is how you add extras into the command for later. You
# can do whatever you like with these though.
my $param_flag = (defined($cf{'params'}{'param_flag'})) ? "--param_flag" : '';
my $param_var = (defined($cf{'params'}{'param_var'})) ? "--param_var ".$cf{'params'}{'param_var'} : '';
if(defined($cf{'params'}{'panic'})){
	die "Oh no! '--param panic' was specified!";
}

# Open up our run file in append mode
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Separate file names into single end and paired end. The array
# gives the file names written to the run file from the previous job
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%cf, @{$cf{'prev_job_files'}});

# Initiate the FastQ encoding type variable.
# Once we've found the encoding for one file, we'll assume all others are the same
my $encoding = 0;

#
# Go through each SINGLE END file and run our command
#
if($se_files && scalar(@$se_files) > 0){
	foreach my $file (@$se_files){

		# Start the clock...
		my $timestart = time;

		# Figure out the encoding if we don't already know
		if(!$encoding){
			($encoding) = CF::Helpers::fastq_encoding_type($file);
		}
		# Format the returned encoding type so that our program can understand it
		my $enc = "";
		if($encoding eq 'phred33' || $encoding eq 'phred64' || $encoding eq 'solexa'){
			$enc = '--'.$encoding;
		}

		# Generate a nice output file name
		my $output_fn = $file."_processed.output";

		# Write our command!
		my $command = "my_command $enc $param_flag $param_var -c $cf{cores} -m $cf{mem} -g $cf{refs}{fasta} -i $file -o $output_fn";
		# Write the command to the log file
		warn "\n###CFCMD $command\n\n";

		# Try to run the command - returns 0 on success (which evaluated to false)
		if(!system ($command)){
			# Command worked!
			# Work out how long the processing took
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);

			# Print a success message to the log file which will be e-mailed out
			warn "###CF Example module (SE mode) successfully exited, took $duration..\n";

			# Check we can find our output filename!
			if(-e $output_fn){
				# Print the current job ID and the output filename to the run file
				# This is so that subsequent modules can use this output
				print RUN "$cf{job_id}\t$output_fn\n";
			} else {
				# Oops - can't find the output file!
				warn "\n###CF Error! Example module output file $output_fn not found..\n";
			}
		} else {
			# Command returned a non-zero result, probably went wrong...
			warn "\n###CF Error! Example module (SE mode) failed for input file '$file': $? $!\n\n";
		}

	}
}


#
# Go through each pair of PAIRED END files and run our command
#
if($pe_files && scalar(@$pe_files) > 0){
	foreach my $files_ref (@$pe_files){
		my @files = @$files_ref;

		# Check that we do actually have two files here
		if(scalar(@files) == 2){

			# Start the clock...
			my $timestart = time;

			# Figure out the encoding if we don't already know
			if(!$encoding){
				($encoding) = CF::Helpers::fastq_encoding_type($files[0]);
			}
			# Format the returned encoding type so that our program can understand it
			my $enc = "";
			if($encoding eq 'phred33' || $encoding eq 'phred64' || $encoding eq 'solexa'){
				$enc = '--'.$encoding.'-quals';
			}

			# Generate a nice output file name
			my $output_fn = $files[0]."_".$files[1]."_processed.output";

			# Write our command!
			my $command = "my_command $enc $param_flag $param_var -c $cf{cores} -m $cf{mem} -g $cf{refs}{fasta} -1 $files[0] -2 $files[1] -o $output_fn";
			# Write the command to the log file
			warn "\n###CFCMD $command\n\n";

			# Try to run the command - returns 0 on success (which evaluated to false)
			if(!system ($command)){
				# Command worked!
				# Work out how long the processing took
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);

				# Print a success message to the log file which will be e-mailed out
				warn "###CF Example module (PE mode) successfully exited, took $duration....\n";

				# Check we can find our output filename
				if(-e $output_fn){
					# Print the current job ID and the output filename to the run file
					# This is so that subsequent modules can use this output
					print RUN "$cf{job_id}\t$output_fn\n";
				} else {
					# Oops - can't find the output file!
					warn "\n###CF Error! Example module output file $output_fn not found\n";
				}
			} else {
				# Command returned a non-zero result, probably went wrong...
				warn "\n###CF Error! Example module (PE mode) failed for input file '$files[0]': $? $!\n\n";
			}
		} else {
			# We didn't have two files here.. This shouldn't ever happen really.
			warn "\n###CF Error! Example module paired end files had ".scalar(@files)." input files instead of 2\n";
		}

	}
}


# Close the run file
close (RUN);

# Any additional cleanup code can come here
# For example, you could move, delete and rename generated files to keep things tidy
