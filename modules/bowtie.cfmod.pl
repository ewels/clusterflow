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

# Save the args for later
my @rawargs = @ARGV;

# Module requirements
my %requirements = (
	'cores' 	=> ['1', '8'],
	'memory' 	=> ['4G', '5G'],
	'modules' 	=> ['bowtie','bowtie2','samtools'],
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# Bowtie alignment typically takes less than 10 hours per BAM file
		# This is probably inaccurate? May need tweaking.
		return CF::Helpers::minutes_to_timestamp ($num_files * 14 * 60);
	}
);

# Help text
my $helptext = "".("-"x22)."\n Bowtie 1 or 2 Module\n".("-"x22)."\n
This module inspects the first input file and calculates the
read length. If this is >= 50bp, it aligns with the Bowtie 2
module. If not, it aligns with the Bowtie 1 module.\n
This module assumes that the bowtie1.cfmod.pl and bowtie2.cfmod.pl
module files are contained within the same directory as this script.\n
See cf --help bowtie1 and cf --help bowtie2 for more information
on these modules.\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);


# MODULE

# Look at the read length of the first file
if(!CF::Helpers::fastq_min_length($cf{'prev_job_files'}[0], 50)){
	warn "\n\n###CF First file has reads < 50bp long. Using bowtie 1 for alignment.\n";

	my $command = "$FindBin::Bin/bowtie1.cfmod.pl ".join(" ", @rawargs);
	warn "\nBowtie 1 module command: $command\n\n";

	system($command);

} else {
	warn "\n\n###CF First file has reads >= 50bp long. Using bowtie 2 for alignment.\n";

	my $command = "$FindBin::Bin/bowtie2.cfmod.pl ".join(" ", @rawargs);
	warn "\nBowtie 2 module command: $command\n\n";

	system($command);

}
