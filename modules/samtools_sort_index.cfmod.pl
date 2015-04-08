#!/usr/bin/env perl
use warnings;
use strict;
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

# Module requirements
my %requirements = (
	'cores' 	=> ['2', '4'],
	'memory' 	=> ['8G', '30G'],
	'modules' 	=> ['samtools'],
	'time' 		=> sub {
		my $runfile = $_[0];
		my $num_files = $runfile->{'num_starting_merged_aligned_files'};
		# Sorting + indexing typically takes less than 1 hour per BAM file
		return $num_files * 60;
	}
);

# Help text
my $helptext = "".("-"x30)."\n Samtools Sort & Index Module\n".("-"x30)."\n
Tries to index BAM / SAM files. If fails, assumes that file is not sorted
and runs samtools sort, then attempts to index again.
Will assume anything not ending in .bam is a SAM file and will convert to
BAM. Output is basename_srtd.bam if sorted.
Index files are written to disk but not written to CF run files.
Using parameter 'byname' or '-n' in pipeline forces sorting by read name.\n";

# Setup
my %runfile = CF::Helpers::module_start(\@ARGV, \%requirements, $helptext);

# MODULE
my $timestart = time;

my $mem = CF::Helpers::human_readable_to_bytes($runfile{'memory'});
my $mem_per_thread = CF::Helpers::bytes_to_human_readable(int($mem / $runfile{'cores'}));
warn "\n\nSamtools memory per thread: $mem_per_thread. Cores: ".$runfile{'cores'}."\n\n\n";

my $namesort = '';
$namesort = '-n' if (defined($runfile{'params'}{'byname'}));

open (RUN,'>>',$runfile{'run_fn'}) or die "###CF Error: Can't write to ".$runfile{'run_fn'}.": $!";

# Print version information about the module.
warn "---------- Samtools version information ----------\n";
warn `samtools 2>&1 | head -n 4`;
warn "------- End of Samtools version information ------\n";

# we want e.g. samtools view -bS ./input.sam | samtools sort - outfile
foreach my $file (@{$runfile{'prev_job_files'}}){

	# Figure out the file type
	my $filetype = "";
	if (lc($file) =~ /\.([sb]am$)/){
		$filetype = $1;
		warn "\n$file looks like a $filetype file\n";
	} else {
		warn "\n Can't determine file-type for $file. Assuming sam... \n";
		$filetype = "sam";
	}

	# Try to index if we have a BAM file
	if($filetype eq "bam"){
		if(samtools_index($file)){
			# samtools worked - print out resulting filenames
			print RUN $runfile{'job_id'}."\t$file\n";
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
			print RUN $runfile{'job_id'}."\t$output_fn\n";

			# Index the sorted file
			if(samtools_index($file)){
				unless (-e "$file.bai"){
					warn "\n###CF Error! samtools index output file $file.bai not found..\n";
				}
			} else {
				warn "\n###CF Error! samtools index failed for $file\n\n";
			}

		} else {
			warn "\n###CF Error! samtools sort output file $output_fn(.bam) not found..\n";
		}

	} else {
		warn "\n###CF Error! samtools sort failed, exited in an error state: $? $!\n\n";
	}
}


my $duration =  CF::Helpers::parse_seconds(time - $timestart);
warn "###CF samtools sort / index finished, took $duration..\n";


# we want e.g. samtools view -bS ./input.sam | samtools sort - outfile
sub samtools_index {
	my $fn = $_[0];
	my $command = "samtools index $fn";
	warn "\n###CFCMD $command\n\n";

	# Return the opposite of the result so that the function returns true if it works
	if(!system ($command)){
		return 1;
	} else {
		return 0;
	}
}



close (RUN);
