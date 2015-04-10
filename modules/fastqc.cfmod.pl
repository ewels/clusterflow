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

# Module requirements
my %requirements = (
	'cores' 	=> '1',
	'memory' 	=> '2G',
	'modules' 	=> ['fastqc'],
	'time' 		=> sub {
		my $runfile = $_[0];
		my $num_files = $runfile->{'num_starting_merged_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# FastQC typically takes less than 12 minutes per file
		return CF::Helpers::minutes_to_timestamp ($num_files * 20);
	}
);

# Help text
my $helptext = "".("-"x15)."\n FastQC Module\n".("-"x15)."\n
FastQC is a quality control tool for high throughput sequence data.
For further information, please run fastqc --help\n\n";

# Setup
my %runfile = CF::Helpers::module_start(\@ARGV, \%requirements, $helptext);


# MODULE
open (RUN,'>>',$runfile{'run_fn'}) or die "###CF Error: Can't write to $runfile{run_fn}: $!";

# Print version information about the module.
warn "---------- FastQC version information ----------\n";
warn `fastqc --version`;
warn "\n------- End of FastQC version information ------\n";

# Read any options from the pipeline parameters
my $nogroup = defined($runfile{'params'}{'nogroup'}) ? "--nogroup" : '';

# Go through each supplied file and run FastQC.
foreach my $file (@{$runfile{'prev_job_files'}}){
	my $timestart = time;
	
	my $command = "fastqc -q $nogroup $file";
	warn "\n###CFCMD $command\n\n";

	if(!system ($command)){
		print RUN "$runfile{job_id}\t$file\n";
		my $duration =  CF::Helpers::parse_seconds(time - $timestart);
		warn "###CF FastQC successfully ran, took $duration\n";
	} else {
		print "###CF Error! FastQC Failed for input file '$file': $? $!\n";
	}
}

close (RUN);
