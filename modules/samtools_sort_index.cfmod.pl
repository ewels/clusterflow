
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
	'modules' 	=> 'samtools',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Sorting + indexing typically takes less than 1 hour per BAM file
		return CF::Helpers::minutes_to_timestamp ($num_files * 60);
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
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

# MODULE
my $mem = CF::Helpers::human_readable_to_bytes($cf{'memory'});
my $mem_per_thread = CF::Helpers::bytes_to_human_readable(int($mem / $cf{'cores'}));
warn "\n\nSamtools memory per thread: $mem_per_thread. Cores: $cf{cores}\n\n\n";

# Set up optional parameters
my $namesort = (defined($cf{'params'}{'byname'})) ? '-n' : '';

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- Samtools version information ----------\n";
warn `samtools 2>&1 | head -n 4`;
warn "------- End of Samtools version information ------\n";

# we want e.g. samtools view -bS ./input.sam | samtools sort - outfile
foreach my $file (@{$cf{'prev_job_files'}}){
	my $timestart = time;

	# Output file name - the .bam is looked for later
	my $output_suffix = "_srtd";
	$output_suffix .= 'n' if ($namesort eq '-n');
	my $output_fn = $file.$output_suffix;

	# Figure out the file type
	my $filetype = "";
	if (lc($file) =~ /\.([sb]am$)/){
		$filetype = $1;
		warn "\n$file looks like a $filetype file\n";
		# Recreate the output filename without the BAM extension.
		($output_fn = $file) =~ s/\.[sb]am$//;
		$output_fn .= $output_suffix;
	} else {
		warn "\n Can't determine file-type for $file. Assuming sam... \n";
		$filetype = "sam";
	}

	# Try to index if we have a BAM file
	warn "Attempting to index input file in case it's already sorted..\n";
	if($filetype eq "bam"){
		# Samtools index returns 0 even if it fails. Look for the .bai file instead.
		my $index_command_one = "samtools index $file 2>&1";
		my $indexing_output  = `$index_command_one`;
		if(-e "$file.bai"){
			warn $indexing_output."\n";
			# samtools worked - print out resulting filenames
			warn "\n###CFCMD $index_command_one\n\n";
			print RUN $cf{'job_id'}."\t$file\n";
			my $duration =  CF::Helpers::parse_seconds(time - $timestart);
			warn "###CF samtools index successfully exited, took $duration. Skipping sort.\n";
			# If we could index, file must already be sorted.
			next;
		} else {
			warn "Samtools index didn't work, file not sorted. Going on to sorting step...\n";
		}
	}

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
			print RUN $cf{'job_id'}."\t$output_fn\n";

			# Index the sorted file
			$timestart = time;
			my $index_command_two = "samtools index $output_fn";
			warn "\n###CFCMD $index_command_two\n\n";
			system ($index_command_two);
			if(-e "$output_fn.bai"){
				my $duration =  CF::Helpers::parse_seconds(time - $timestart);
				warn "###CF samtools index successfully exited, took $duration.\n";
			} else {
				warn "\n###CF Error! samtools index failed for $file\n\n";
			}

		} else {
			warn "\n###CF Error! samtools sort output file $output_fn not found..\n";
		}

	} else {
		warn "\n###CF Error! samtools sort failed, exited in an error state: $? $!\n\n";
	}
}


close (RUN);
