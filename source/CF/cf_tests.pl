#!/usr/bin/env perl
use warnings;
use strict;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptionsFromArray);
use IPC::Open3;
use lib "$FindBin::Bin/../";
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

# Variable setup
my $cf_version = $CF::Constants::CF_VERSION;
my $homedir = $ENV{"HOME"};
my %config = %CF::Constants::config;
my $base_dir = "$Bin/../../";
my @pipeline_folders = ("$homedir/clusterflow/pipelines/", "$Bin/../../pipelines/");
my @module_folders = ("$homedir/clusterflow/modules/", "$Bin/../../modules/");

# Command line options
my $version;
my $help;
my $config_result = GetOptions(
	"version" => \$version,
	"help" => \$help
);

# Info command line calls
if($version){
	print "\nCluster Flow v$cf_version\n\n";
	exit;
}
if($help){
    print "".("-"x25)."\n Cluster Flow Test Suite\n".("-"x25)."\nCluster Flow v$cf_version\n\n";
    print "This is a stand-alone script with a test-suite for Cluster Flow.\n".
        "It is designed to aid development by making it easier to check\n".
        "for bugs in the code before pushing updates.\n\n";
    exit;
}


#
# STEP ONE - CORE CF SYNTAX
#
print "## Step One: Checking core Cluster Flow syntax.\n";
system("perl -c ".$base_dir."cf");


#
# STEP TWO - MODULE SYNTAX CHECKS
#
print "\n## Step Two: Checking module script syntax.\n";
my $num_passed = 0;
my $num_failed = 0;
my $unrecognised_filetype = 0;
my @modules;
foreach my $folder (@module_folders){
    if(-e $folder){
        opendir (DIR, $folder) or die $!;
        my @dir_files = sort readdir(DIR);
        while ( my $file = shift @dir_files ) {
            my $path = $folder.$file;
            # Check Perl files
            if($file =~ /\.pl$/){
                # Crude way to hide the output from valid files
                if(system("perl -c $path > /dev/null 2>&1")){
                    system("perl -c $path");
                    $num_failed++;
                } else {
                    push(@modules, $path);
                    $num_passed++;
                }
            }
            # Check Python files
            elsif($file =~ /\.py$/){
                if(system("python -m py_compile $path > /dev/null 2>&1")){
                    print "Error! $path did not compile\n";
                    $num_failed++;
                } else {
                    push(@modules, $path);
                    $num_passed++;
                }
            }
            # Unreocognised file type
            elsif($file ne '.' && $file ne '..') {
                $unrecognised_filetype++;
            }
        }
        closedir(DIR);
    }
}
print "$num_passed modules passed, $num_failed failed and $unrecognised_filetype had an unrecognised file type.\n";

#
# STEP THREE - TRY TO GET REQUIREMENTS FROM PASSED MODULES
#
print "\n## Step Three: Get requirements from $num_passed passed modules.\n";
my $runfn = "$Bin/../../clusterflow.config"; # cheatin'
my $num_requirements_failed = 0;
foreach my $module_fn (@modules){

    # Stolen and modified from core cf - note, needs to be kept up to date
    my $failed = 0;

	# Send query to module
	my $response = `$module_fn --requirements --run_fn $runfn --cores 6 --mem 128G 2> /dev/null`;

	# Parse response
	my $required_modules;
	my $cores = 'unset';
	my $mem = 'unset';
	my $time = 'unset';
	foreach my $ln (split("\n", $response)){
		chomp($ln);
		my ($key, $val) = split(':', $ln, 2);
		if($key and $val){
			# trim whitespace
			$key =~ s/^\s+|\s+$//g;
			$val =~ s/^\s+|\s+$//g;

			if($key eq 'cores'){
				$cores = $val;
				$cores =~ s/\D//g;
			}
			if($key eq 'memory'){
				$mem = $val;
			}
			if($key eq 'modules'){
				my @modules = split(/[\s,]+/, $val);
			}
			if($key eq 'time'){
				$time = $val;
			}
		}
	}

	# Check that we have what we need
    if($cores eq 'unset'){
        print "$module_fn cores not set\n";
        $failed = 1;
    } else {
        unless($cores =~ /^\d+$/){
            print "$module_fn cores not numeric ('$cores')\n";
            $failed = 1;
        }
    	if($cores > 6){
            print "$module_fn cores greater than 6 ('$cores') \n";
            $failed = 1;
        }
    }


    if($mem eq 'unset'){
        print "$module_fn memory not set\n";
        $failed = 1;
    } else {
        unless($mem =~ /^[\d\.]+[gmkb]*$/i){
            print "$module_fn memory wrong format ('$mem') \n";
            $failed = 1;
        }
    	if(CF::Helpers::human_readable_to_bytes($mem) > CF::Helpers::human_readable_to_bytes('128G')){
            print "$module_fn memory > allocated 128G ('$mem')\n";
            $failed = 1;
        }
    }


    if($time eq 'unset'){
        print "$module_fn time not set\n";
        $failed = 1;
    } else {
        unless(CF::Helpers::timestamp_to_minutes($time) > 0){
            print "$module_fn time not recognised ('$time')\n";
            $failed = 1;
        }
    }

    if($failed){
        $num_requirements_failed++;
        print "\nCommand for testing:\n$module_fn --requirements --run_fn $runfn --cores 6 --mem 128G 2> /dev/null\n";
        print "".('-'x50)."\n";
    }
}
print "".($num_passed - $num_requirements_failed)." modules passed requirements check, $num_requirements_failed failed.\n";

print "\nDone.\n\n";
