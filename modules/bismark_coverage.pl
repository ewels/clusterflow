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
	'memory' 	=> ['3G', '10G'],
	'modules' 	=> 'bismark',
	'time' 		=> sub {
		my $runfile = $_[0];
		my $num_files = $runfile->{'num_starting_merged_aligned_files'};
        $num_files = ($num_files > 0) ? $num_files : 1;
		# Bismark coverage analysis typically takes less than 2 hours per BAM file
		return CF::Helpers::minutes_to_timestamp ($num_files * 4 * 60);
	}
);

# Help text
my $helptext = "".("-"x38)."\n Bismark Coverage Analysis Module\n".("-"x38)."\n
The bismark_coverage module runs the bismark coverage2cytosine script.
This reads in a fasta reference and creates a .cov file with coverage information
about every cytosine in the genome. The module then uses this to plot two
informative plots about the coverage of the experiment.
Use the capture_regions=<filename.txt> paramter when doing a sequence capture
experiment for information about targeted enrichment.\n".
"Use coverage2cytosine --help for more information.\n\n";

# Setup
my %runfile = CF::Helpers::module_start(\@ARGV, \%requirements, $helptext);

# MODULE
my $timestart = time;

# Check that we have a genome defined
if(!defined($runfile{'config'}{'references'}{'fasta'})){
   die "\n\n###CF Error: No genome fasta path found in run file $runfile{run_fn} for job $runfile{job_id}. Exiting..";
} else {
    warn "\nAligning against $runfile{config}{references}{fasta}\n\n";
}

open (RUN,'>>',$runfile{'run_fn'}) or die "###CF Error: Can't write to $runfile{run_fn}: $!";

# Print version information about the module.
warn "---------- bismark coverage2cytosine version information ----------\n";
warn `coverage2cytosine --version`;
warn "\n------- End of bismark_methylation_extractor version information ------\n";

# Read any options from the pipeline parameters
my $capture_regions = (defined($runfile{'params'}{'capture_regions'})) ? $runfile{'params'}{'capture_regions'} : 0;


# Go through each file and deduplicate
foreach my $file (@{$runfile{'prev_job_files'}}){

    #
    # Make the genome wide coverage report
    #
	# Should have been given a deduplicated BAM file
	$file =~ s/.[bs]am$//;
	$file .= '.bismark.cov';
    my $output_cov_fn = substr($file,0 ,-3)."gwCov.cov";
    my $cmd = "coverage2cytosine --genome_folder $runfile{config}{references}{fasta} $file -o $output_cov_fn";
    warn "\n###CFCMD $cmd\n\n";

    if(!system ($cmd)){
        # Bismark worked - print out resulting filenames
        my $duration =  CF::Helpers::parse_seconds(time - $timestart);
        warn "###CF Bismark coverage2cytosine successfully exited, took $duration..\n";
        if(-e $output_cov_fn){
            print RUN "$runfile{job_id}\t$output_cov_fn\n";
        } else {
            warn "\n###CF Error! Bismark coverage2cytosine output file $output_cov_fn not found..\n";
            next;
        }
    } else {
        warn "\n###CF Error!Bismark coverage2cytosine exited with an error state for file '$file': $? $!\n\n";
    }


	# perl script paths
	my $coverage_script = "/home/phil/scripts/ngi_visualizations/stand_alone/bismark/bismark_coverage_curves.pl";
	my $windows_script = "/home/phil/scripts/ngi_visualizations/stand_alone/bismark/bismark_window_sizes.pl";

	my $regions = '';
	if($capture_regions and -e $capture_regions){
		$regions = "--regions $capture_regions";
	} elsif($capture_regions){
		warn "###CF Error! Can't find capture regions BED file '$capture_regions' - running without capture information.\n";
	}

	# Plot coverage curves
	if(-e $coverage_script){
		my $cov_cmd = "perl $coverage_script $regions --stranded $output_cov_fn";
		warn "###CFCMD $cov_cmd\n\n";
		if(!system ($cov_cmd)){
			warn "###CF Coverage curves plot successfully created.\n";
		} else {
			warn "###CF Error! Coverage curves plotting exited with an error code for input file '$output_cov_fn'.\n";
		}
	} else {
		warn "###CF Warning: Couldn't find find the coverage script: $coverage_script\nSkipping..\n";
	}
	# Plot coverage curves
	if(-e $windows_script){
		my $win_cmd = "perl $windows_script $regions $output_cov_fn";
		warn "###CFCMD $win_cmd\n\n";
		if(!system ($win_cmd)){
			warn "###CF Window sizes plot successfully created.\n";
		} else {
			warn "###CF Error! Window size plotting exited with an error code for input file '$output_cov_fn'.\n";
		}
    } else {
		warn "###CF Warning: Couldn't find find the windows script: $windows_script\nSkipping..\n";
	}
}
