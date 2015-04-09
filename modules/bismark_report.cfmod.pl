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
	'memory' 	=> '100M',
	'modules' 	=> 'bismark',
	'time' 		=> '10'
);

# Help text
my $helptext = "".("-"x21)."\n Bismark Report Module\n".("-"x21)."\n
This script runs the bismark2report script to generate an
overview report. It will run on everything in the directory,
overwriting previously generated reports.\n
For further information, please run bismark2report --help\n\n";

# Setup
my %runfile = CF::Helpers::module_start(\@ARGV, \%requirements, $helptext);


# --version. Returns version information about the module.
warn "---------- bismark2report version information ----------\n";
warn `bismark2report --version`;
warn "\n------- End of bismark2report version information ------\n";

if(!system ("bismark2report")){
	warn "###CF Bismark report successfully created\n";
} else {
	warn "###CF Error! Bismark report exited with an error state: $? $!\n";
}
