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
#        the filename to MODULENAME.cfmod(.pl)
#

#
# Cluster Flow Job Initialisation
#

# Get options
my $avail_cores;
my $avail_mem;
my $required_modules;
my $run_fn;
my $help;
my $result = GetOptions (
	"cores=i" => \$avail_cores,
	"mem=s" => \$avail_mem,
	"modules" => \$required_modules,
    "runfn" => \$run_fn,
	"help" => \$help
);

# --cores specifies how many cores are recommended
# Return how many cores we need. This can be more or less than
# the number recommended.
if($avail_cores){
	# The allocate_cores() function looks at what's been offered,
	# compares this to a minimum and maximum (1 and 6 in this case)
	# and returns a valid number
	print CF::Helpers::allocate_cores($avail_cores, 1, 6);
	exit;
}

# --mem specifies how much memory is recommended
# Return how much we need. This can be more or less than the amount recommended.
if($avail_mem){
	# The allocate_memory() function looks at what's been offered,
	# compares this to a minimum and maximum (3G and 4G in this case)
	# and returns a valid number. Works with any human readable string.

	my ($runfn, $cores, @parameters) = @ARGV;
	my ($starting_files, $config) = CF::Helpers::parse_runfile_prerun($runfn);
	
	# examples of how to extract extra information to make memory decision
	# warn "cores: $cores";
	# warn "param: $_" foreach(@parameters);  # display params for this line of pipeline
	# warn "config key: $_  value: $$config{$_}" foreach(keys(%$config));
	# warn "starting file: $_" foreach(@$starting_files);
	# warn "available mem: $avail_mem\n";
	
	# return memory request
	print CF::Helpers::allocate_memory($avail_mem, '3G', '4G');	
	exit;
}

# --modules. Return comma seperated names of any environment modules which should be loaded.
if($required_modules){
	print 'example_module,sratoolkit';
	exit;
}

# --help. Print the help.
if($help){
	print "".("-"x15)."\n Example Perl Module\n".("-"x15)."\n
This is a dummy module which provides example code to be used
by anyone interested in writing their own Cluster Flow module.
It doesn't do anything, but its source code is nicely commented..\n\n";
	exit;
}



#
# PROPER MODULE CODE BELOW
#


# Start the clock...
my $timestart = time;

# Read in the input files from the run file
my (
	$files,
	$runfile,
	$job_id,
	$prev_job_id,
	$cores,
	$mem,
	$parameters,
	$config_ref
) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;

# Check that we have a reference genome defined, die if we don't.
if(!defined($config{references}{fasta})){
   die "\n\n###CF Error: No genome fasta path found in run file $runfile for job $job_id. Exiting.. ###";
} else {
    warn "\nUsing the genome path: ".$config{references}{fasta}."\n\n";
}

# Sanity check for cores and memory, shouldn't be required but can't hurt
if(!defined($cores) || $cores < 1){
	$cores = 1;
}
if(CF::Helpers::human_readable_to_bytes($mem) < CF::Helpers::human_readable_to_bytes('3G')){
	$mem = '3G';
}

# Print version information about the program to be executed.
warn "---------- < module > version information ----------\n";
warn `MY_COMMAND --version`;
warn "\n------- End of < module > version information ------\n";	

# Initiate the FastQ encoding type variable.
# Once we've found the encoding for one file, we'll assume all others are the same
my $encoding = 0;

# Read any options from the pipeline parameters
# Prepare empty variables that can contain extra command line options
my $param1 = "";
my $param2 = "";
foreach my $parameter (@$parameters){
	if($parameter eq "option_1"){
		$param1 = "--parameter1";
	}
	if($parameter eq "option_2"){
		$param2 = "--parameter2";
	}
	if($parameter eq "option_3"){
		$param1 = "--parameter1";
		$param2 = "--parameter2";
	}
}

# Open up our run file in append mode
open (RUN,'>>',$runfile) or die "###CF Error: Can't write to $runfile: $!";

# Separate file names into single end and paired end
my ($se_files, $pe_files) = CF::Helpers::is_paired_end(\%config, @$files);

#
# Go through each SINGLE END file and run our command
#
if($se_files && scalar(@$se_files) > 0){
	foreach my $file (@$se_files){
		
		# Figure out the encoding if we don't already know
		if(!$encoding){
			($encoding) = CF::Helpers::fastq_encoding_type($file);
		}
		# Format the returned encoding type so that our program can understand it
		my $enc = "";
		if($encoding eq 'phred33' || $encoding eq 'phred64' || $encoding eq 'solexa'){
			$enc = '--'.$encoding.'-quals';
		}
		
		# Generate a nice output file name
		my $output_fn = $file."_processed.output";
		
		# Write our command!
		my $command = "my_command $enc $param1 $param2 -c $cores -m $mem -g ".$config{references}{fasta}." -i $file -o $output_fn";
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
				print RUN "$job_id\t$output_fn\n"; 
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
			my $command = "my_command $enc $param1 $param2 -c $cores -m $mem -g ".$config{references}{fasta}." -1 ".$files[0]." -2 ".$files[1]." -o $output_fn";
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
					print RUN "$job_id\t$output_fn\n";
				} else {
					# Oops - can't find the output file!
					warn "\n###CF Error! Example module output file $output_fn not found\n";
				}
			} else {
				# Command returned a non-zero result, probably went wrong...
				warn "\n###CF Error! Example module (PE mode) failed for input file '".$files[0]."': $? $!\n\n";
			}
		} else {
			# We didn't have two files here..
			warn "\n###CF Error! Example module paired end files had ".scalar(@files)." input files instead of 2\n";
		}
		
	}
}


# Close the run file
close (RUN);

# Any additional cleanup code can come here
# For example, you could move, delete and rename generated files to keep things tidy




