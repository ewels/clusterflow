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
	'modules' 	=> ['multiqc'],
	'time' 		=> '30'
);

# Help text
my $helptext = "".("-"x15)."\n MultiQC\n".("-"x15)."\n
MultiQC creates summary reports showing analysis results
across multiple samples for any analysis pipeline.
For further information, please run multiqc --help\n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);


# MODULE
open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- MultiQC version information ----------\n";
warn `multiqc --version`;
warn "\n------- End of MultiQC version information ------\n";

# Read any options from the pipeline parameters
my $template = defined($cf{'params'}{'template'}) ? "-t ".$cf{'params'}{'template'} : '';

# Run MultiQC - doesn't need to know about previous files
my $timestart = time;

my $command = "multiqc -v $template .";
warn "\n###CFCMD $command\n\n";

if(!system ($command)){
	print RUN "$cf{job_id}\t$file\n";
	my $duration =  CF::Helpers::parse_seconds(time - $timestart);
	warn "###CF MultiQC successfully ran, took $duration\n";
} else {
	print "###CF Error! MultiQC exited with an error state: $? $!\n";
}

close (RUN);
