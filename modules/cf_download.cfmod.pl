#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;
use POSIX;

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
	'memory' 	=> '1G',
	'modules' 	=> '',
	'time' 		=> sub {
		my $runfile = $_[0];
		my $num_files = $runfile->{'num_starting_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
		# This is a tough one! Let's assume a dial-up connection..
		return CF::Helpers::minutes_to_timestamp ($num_files * 4 * 60);
	}
);

# Help text
my $helptext = "\nThis is a core module which downloads files for Cluster Flow.\n";

# Setup
my %runfile = CF::Helpers::module_start(\@ARGV, \%requirements, $helptext);

# MODULE
my $timestart = time;

# Strip number from download job ID so that these files are all read in the next module
my $job_id = $runfile{'job_id'};
$job_id =~ s/download_[\d]{3}$/download/;

# Get URL and download filename from params
my $url = $runfile{'params'}{'url'};
my $dl_fn = $runfile{'params'}{'dl_fn'};

warn "\n---------------------\nDownloading $dl_fn from $url\nStarted at ".strftime("%H:%M, %A - %d/%m/%Y", localtime)."\n";

open (RUN,'>>',$runfile{'run_fn'}) or die "###CF Error: Can't write to $runfile{run_fn}: $!";

my $command = "wget -nv --tries=10 --output-document=$dl_fn $url";
warn "\n###CFCMD $command\n\n";
if(!system ($command)){
	# Download worked - print resulting filename to results file
	print RUN "$job_id\t$dl_fn\n";
	my $duration =  CF::Helpers::parse_seconds(time - $timestart);
	warn "###CF Download worked - took $duration\n";
} else {
	# Download failed - don't print a filename so that child processes exit silently
	warn "###CF Error! Download '$dl_fn' failed: $? $!\n";
}

my $date = strftime "%H:%M %d-%m-%Y", localtime;
warn "\nDownload module finished at $date\n---------------------\n";

close (RUN);
