#!/usr/bin/env perl
use warnings;
use strict;
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../source";
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
	'cores' 	=> '1',
	'memory' 	=> '2G',
	'modules' 	=> 'samtools',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Converting typically takes less than 10 minutes per BAM file
		return CF::Helpers::minutes_to_timestamp ($num_files * 20);
	}
);

# Help text
my $helptext = "".("-"x30)."\n Samtools bam2sam Module\n".("-"x30)."\n
Converts BAM to SAM files using samtools view.\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);

open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- Samtools version information ----------\n";
warn `samtools 2>&1 | head -n 4`;
warn "------- End of Samtools version information ------\n";

# we want e.g. samtools view -bS ./input.sam | samtools sort - outfile
foreach my $file (@{$cf{'prev_job_files'}}){
	my $timestart = time;

	# Output file name
	(my $output_fn = $file) =~ s/\.bam$//;
	$output_fn .= '.sam';

	# Command
	my $cmd .= "samtools view $file > $output_fn";
	warn "\n###CFCMD $cmd\n\n";

	if(!system ($cmd)){
		# samtools worked - print out resulting filenames
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF samtools view successfully exited, took $duration..\n";

		if(-e $output_fn){
			print RUN $cf{'job_id'}."\t$output_fn\n";
		} else {
			warn "\n###CF Error! samtools view output file $output_fn not found..\n";
		}

	} else {
		warn "\n###CF Error! samtools view failed, exited in an error state: $? $!\n\n";
	}
}


close (RUN);
