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
my $result = GetOptions ("cores=i" => \$required_cores, "mem=s" => \$required_mem, "modules" => \$required_modules, "runfn=s" => \$run_fn, "help" => \$help);

# QSUB SETUP
# --cores i = offered cores. Return number of required cores.
if($required_cores){
	print CF::Helpers::allocate_cores($required_cores, 2, 4);
	exit;
}
# --mem. Return the required memory allocation.
if($required_mem){
    print CF::Helpers::allocate_memory($required_mem, '8G', '30G');
}
# --modules. Return csv names of any modules which should be loaded.
if($required_modules){
	print 'samtools';
	exit;
}
# --help. Print help.
if($help){
	print "".("-"x23)."\n Samtools Sort & Index Module\n".("-"x23)."\n
Tries to index BAM / SAM files. If fails, assumes that file is not sorted
and runs samtools sort, then attempts to index again.
Will assume anything not ending in .bam is a SAM file and will convert to
BAM. Output is basename_srtd.bam if sorted.
Index files are written to disk but not written to CF run files.
Using parameter 'byname' or '-n' in pipeline forces sorting by read name.\n";
	exit;
}

# MODULE
my $timestart = time;

# Read in the input files from the run file
my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config_ref) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;
$mem = CF::Helpers::human_readable_to_bytes($mem);
my $mem_per_thread = int($mem/$cores);
warn "\n\n Samtools: mem per thread: $mem_per_thread ; cores: $cores\n\n\n";

if(!defined($cores) or $cores < 1){
	$cores = 1;
}

my $namesort = '';
$namesort = '-n' if (grep(/^byname$/, @$parameters) > 0);

open (RUN,'>>',$runfile) or die "###CF Error: Can't write to $runfile: $!";

# Print version information about the module.
warn "---------- Samtools version information ----------\n";
warn `samtools 2>&1 | head -n 4`;
warn "\n------- End of Samtools version information ------\n";

# we want e.g. samtools view -bS ./input.sam | samtools sort - outfile
if($files && scalar(@$files) > 0){
	foreach my $file (@$files){

		# Figure out the file type
		my $filetype = "";
		if (lc($file) =~ /\.([sb]am$)/){
			$filetype = $1;
			warn "\nGuessing file $file is a $filetype file\n";
		} else {
			warn "\n Can't determine file-type for $file. Assuming sam... \n";
			$filetype = "sam";
		}

		# Try to index if we have a BAM file
		if($filetype eq "bam"){
			if(samtools_index($file)){
				# samtools worked - print out resulting filenames
				print RUN "$job_id\t$file\n";
				unless (-e "$file.bai"){
					warn "\n###CF Error! samtools index output file $file.bai not found..\n";
				}
				# If we could index, file must already be sorted.
				next;
			}
		}

		# Output file name
		my $output_fn = $file."_srtd";
		$output_fn .= 'n' if ($namesort eq '-n');

		# Pipe BAM stream if we need it
		my $command = '';
		my $sortfile = $file;
		if ($filetype eq "sam"){
			$command .= "samtools view -bS -u $file | ";
			$sortfile = "-";
		}

		$command .= "samtools sort -m $mem_per_thread $namesort $sortfile $output_fn";
		warn "\n###CFCMD $command\n\n";

		if(!system ($command)){
			# samtools worked - print out resulting filenames
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF samtools sort successfully exited, took $duration..\n";

			# Did we add .bam to the filename?
			if (-e "$output_fn.bam"){
				$output_fn = "$output_fn.bam";
			}

			if(-e $output_fn){
				print RUN "$job_id\t$output_fn\n";

				# Index the sorted file
				if(samtools_index($file)){
					unless (-e "$file.bai"){
						warn "\n###CF Error! samtools index output file $file.bai not found..\n";
					}
				} else {
					warn "\n###CF Error! samtools index failed for $file, exited in an error state: $? $!\n\n";
				}

			} else {
				warn "\n###CF Error! samtools sort output file $output_fn(.bam) not found..\n";
			}

		} else {
			warn "\n###CF Error! samtools sort failed, exited in an error state: $? $!\n\n";
		}
	}
}


my $duration =  CF::Helpers::parse_seconds(time - $timestart);
warn "###CF samtools sort / index successfully exited, took $duration..\n";


# we want e.g. samtools view -bS ./input.sam | samtools sort - outfile
sub samtools_index {
	my $command = "samtools index $_";
	warn "\n###CFCMD $command\n\n";

	# Return the opposite of the result so that the function returns true if it works
	if(!system ($command)){
		return 1;
	} else {
		return 0;
	}
}



close (RUN);
