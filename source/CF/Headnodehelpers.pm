#!/usr/bin/perl
package CF::Headnodehelpers;

use warnings;
use strict;
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/source";
use Exporter;
use POSIX qw(strftime);
use XML::Simple;
use Time::Local;
use Term::ANSIColor;
use CF::Helpers;
use Data::Dumper;
use LWP::Simple qw($ua get);

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

# Set the default timeout for LWP::Simple
$ua->timeout(5);

# Function to parse squeue results on a SLURM system
sub parse_localjobs {

	my ($all_users, $cols) = @_;

    # Get job list
    my $osx = 0;
    my $jobs = `ps -o pid -o cmd 2>&1`;
    chomp($jobs);

    # Mac, not linux
    if($jobs =~ /cmd: keyword not found/){
        $osx = 1;
        $jobs = `ps -o pid -o command`;
        chomp($jobs);
    }

	my $output = "";
	foreach my $job (split(/\n/, $jobs)){
        # Look for root cluster flow bash jobs
        if($job =~ /^(\d+) bash cf_local_cf_(.+)_(\d{10})_commands.sh$/){
            # Get parts
            my $pid = $1;
            my $pipeline = $2;
            my $pipelinekey = "$2_$3";
            my $started = $3;

    		# Figure out the cwd for this pipeline
            my $pipeline_wd = `lsof -a -p $pid -d cwd -Fn | tail -1 | sed 's/.//'`;
            $pipeline_wd =~ s/\n//;

    		$output .= "\n".('=' x 70)."\n";
    		$output .= sprintf("%-24s%-44s\n", " Cluster Flow Pipeline:", $pipeline);
    		$output .= sprintf("%-24s%-44s\n", " Top-Level Process ID:", $pid);
    		$output .= sprintf("%-24s%-44s\n", " Submitted:", CF::Helpers::parse_seconds(time - $started)." ago");
    		$output .= sprintf("%-24s%-44s\n", " Working Directory:", $pipeline_wd) if length($pipeline_wd) > 0;
    		$output .= sprintf("%-24s%-44s\n", " Cluster Flow ID:", $pipelinekey);
    		$output .= "".('=' x 70)."\n";

            # Look at child processes using pgrep
            # Should only be one module running
            my $curr_cmd = '?';
            my $curr_mod = '?';
            foreach my $child_pid (split(/\n/, `pgrep -P $pid`)){
                my $child_cmd = '';
                if($osx){
                    $child_cmd = `ps -o command $child_pid`;
                } else {
                    $child_cmd = `ps -o cmd $child_pid`;
                }
                $child_cmd =~ s/[\n\r]//;
                $child_cmd =~ s/COMMAND//;
                $child_cmd =~ s/CMD//;
                $child_cmd =~ s/perl //;
                if($child_cmd =~ /\/modules\/(.+)\.cfmod/){
                    $curr_cmd = $child_cmd;
                    $curr_mod = $1;
                }
            }

            # Look into the bash script to see how many steps there are
            # Count how far through we are.
            if($pipeline_wd){
                my $script_fn = "$pipeline_wd/cf_local_cf_${pipelinekey}_commands.sh";
                my $num_commands = 0;
                my $pos = 0;
                my @mod_names;
                if(-e $script_fn){
                	open (BASH, $script_fn) or die "Can't read $script_fn: $!";
                	while (<BASH>) {
                		chomp;
                		s/[\n\r]//;
                        s/ 2>&1 \| sed .+$//;
                        if(/\/modules\/(.+)\.cfmod/){
                            next if(/cf_run_finished/);
                            next if(/cf_runs_all_finished/);
                            $num_commands++;
                            push (@mod_names, $1);

                            my $this = $_;
                            my $match = $curr_cmd;
                            $this =~ s/\s//g;
                            $match =~ s/\s//g;
                            if($this eq $match){
                                $pos = $num_commands;
                            }
                        }
                    }
                    close(BASH);
                    $output .= "\nCurrently running job $pos out of $num_commands:\n";
                    my $i = 0;
                    foreach my $mod (@mod_names){
                        $i++;
                        if ($i < $pos){
                            if($cols){
                                $output .= color 'green';
                                $output .= "  $mod ";
                                $output .= color 'reset';
                                $output .= "\n"
                            } else {
                                $output .= "  $mod (done)\n";
                            }
                        }
                        elsif($i == $pos){
                            if($cols){
                                $output .= color 'yellow';
                                $output .= "  $mod ";
                                $output .= color 'reset';
                                $output .= "\n"
                            } else {
                                $output .= "> $mod\n";
                            }
                        } elsif ($i > $pos){
                            $output .= "  $mod\n";
                        }
                    }
                    $output .= "\n\n";
                }
            }
        }
	}
	return ($output);
}

# Function to parse squeue results on a SLURM system
sub parse_squeue {

	my ($all_users, $cols) = @_;

	# Set up variables
	my %jobs;
	my %pipelines;
	my @jlist;
	my %pipeline_single_job_ids;

	# Build command
	chomp(my $curruser = `whoami`);
	my $squeue_command = 'squeue -o "%i %T %u %C %S %Q %j %E %r"';
	unless($all_users){
		$squeue_command .= " -u $curruser";
	}
	#TODO - change code so that ordering isn't important
	# For now - sort by job ID
	$squeue_command .= ' |  sort';

	# Parse results
	my $squeue = `$squeue_command`;
	my @sjobs = split(/[\n\r]+/, $squeue);
	foreach my $job (@sjobs){
		# Set up vars for this job
		my %jobhash;
		my $pipeline = 'unknown';
                my $pipelinekey = 'unknown';

		# Split result
		my ($jid, $state, $owner, $cores, $started, $priority, $jobname, $dependency, $dependency_reason) = split(/ /, $job);
		next if($jid eq 'JOBID');

		# Parse
		if($jobname =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
			$pipeline = $1;
			$pipelinekey = "$1_$2";
			$pipelines{$pipelinekey}{started} = $2;
			$jobhash{module} = $3;
		}

		$pipeline_single_job_ids{$pipelinekey} = $jid;

		$jobhash{pipeline} = $pipeline;
		$jobhash{pipelinekey} = $pipelinekey;
		$jobhash{jobname} = $jobname;
		$jobhash{full_job_name} = $jobname;
		$jobhash{state} =  $state;
		$jobhash{cores} = $cores;
		# $jobhash{mem} = ??
		$jobhash{owner} = $owner;
		$jobhash{priority} = $priority;
		$jobhash{started} = $started;
		$jobhash{dependency_reason} = $dependency_reason;
		$jobhash{children} = {};
		$jobhash{parent} = '';

		# If we're queued with dependencies, parse them
		if($dependency_reason eq 'Dependency'){
			# remove any non-numeric stuff, eg. afterany:
			$dependency =~ s/[^\d,]+//g;
			# Split by commas
			my @parents = split(',', $dependency);
			# If more than one, take job with highest ID
			@parents = sort { $a <=> $b } @parents;
			my $parent = $parents[-1];
			# Parse tree
			if($parent && length($parent) > 0){
				parse_job_dependency_hash(\%jobs, $parent, $jid, \%jobhash);
			}
		# No dependencies - give to top level
		} else {
			$jobs{$jid} = \%jobhash;
		}

	}

	# Work out pipeline cwds
	# Figure out the cwd for this pipeline
	my %pipeline_wds;
	foreach my $pipelinekey (keys (%pipelines)){
		my $pipeline_wd;
		if($pipeline_single_job_ids{$pipelinekey}){
			my $cwd_command = "scontrol show job $pipeline_single_job_ids{$pipelinekey} | grep WorkDir";
			$pipeline_wd = `$cwd_command`;
			$pipeline_wd =~ s/\s*WorkDir=//;
			$pipeline_wd =~ s/\n//;
			if(length($pipeline_wd) > 0){
				$pipeline_wds{$pipelinekey} = $pipeline_wd;
			}
		}
	}

	#print Dumper (\%jobs); exit;
	my $output = print_jobs_output(\%jobs, \%pipelines, \%pipeline_single_job_ids, \%pipeline_wds, $cols, $all_users);


	return ($output);

}

# Function to parse qstat results and return them in a nicely formatted manner
sub parse_qstat {

	my ($all_users, $cols) = @_;

	my $qstat_command = "qstat -pri -r -xml";
	if($all_users){
		$qstat_command .= ' -u "*"';
	}

	if (eval {require XML::Parser}) {
			$ENV{XML_SIMPLE_PREFERRED_PARSER} = 'XML::Parser';
	}
	else {
			warn "XML::Parser not available - cf --qstat could be slow\n";
	}

	my $qstat = `$qstat_command`;
	my $xml = new XML::Simple;
	my $data = $xml->XMLin($qstat);

#	print Dumper $data; exit;

	my %jobs;
	my %pipelines;
	my @jlist;
	my %pipeline_single_job_ids;

	# Running Jobs
	# Returns a hash instead of an array if only one element
	if(ref($data->{queue_info}->{job_list}) eq 'HASH'){
		@jlist = \%{$data->{queue_info}->{job_list}};
	} elsif(ref($data->{queue_info}->{job_list}) eq 'ARRAY'){
		@jlist = @{$data->{queue_info}->{job_list}};
	}
	if($data->{queue_info}->{job_list}){
		foreach my $job (@jlist){

			my $jid = $job->{JB_job_number};
			my $pipeline = 'unknown';
			my $pipelinekey = 'unknown';
			my $jobname = $job->{full_job_name};

			if($jobname =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
				$pipeline = $1;
				$pipelinekey = "$1_$2";
				$pipelines{$pipelinekey}{started} = $2;
				$jobs{$jid}{module} = $3;
			}

			$pipeline_single_job_ids{$pipelinekey} = $jid;

			$jobs{$jid}{pipeline} = $pipeline;
			$jobs{$jid}{pipelinekey} = $pipelinekey;
			$jobs{$jid}{jobname} = $jobname;
			$jobs{$jid}{full_job_name} = $job->{full_job_name};
			$jobs{$jid}{state} =  $job->{state}->[0];
			if($job->{state}->[1] eq 'dr' || $job->{state}->[1] eq 't'){
				$jobs{$jid}{state} = 'deleting';
			}
			$jobs{$jid}{cores} = $job->{slots};
			$jobs{$jid}{mem} = $job->{hard_request}->{content};
			$jobs{$jid}{owner} = $job->{JB_owner};
			$jobs{$jid}{priority} = $job->{JB_priority};
			$jobs{$jid}{started} = $job->{JAT_start_time};
			$jobs{$jid}{dependency_reason} = '';
			$jobs{$jid}{children} = {};

		}
	}

	# Pending Jobs
	# Returns a hash instead of an array if only one element
	@jlist = ();
	if(ref($data->{job_info}->{job_list}) eq 'HASH'){
		@jlist = \%{$data->{job_info}->{job_list}};
	} elsif(ref($data->{job_info}->{job_list}) eq 'ARRAY'){
		@jlist = @{$data->{job_info}->{job_list}};
	}
	if($data->{job_info}->{job_list}){
		foreach my $job (@jlist){
			my $jid = $job->{JB_job_number};
			my %jobhash;
			my $pipeline = 'unknown_pending';
			my $pipelinekey = 'unknown_pending';
			my $jobname = $job->{full_job_name};

			if($jobname =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
				$pipeline = $1;
				$pipelinekey = "$1_$2";
				$pipelines{$pipelinekey}{started} = $2;
				$jobhash{module} = $3;
			}

			$pipeline_single_job_ids{$pipelinekey} = $jid;

			$jobhash{pipeline} = $pipeline;
			$jobhash{pipelinekey} = $pipelinekey;
			$jobhash{jobname} = $jobname;
			$jobhash{full_job_name} = $job->{full_job_name};
			$jobhash{state} = $job->{state}->[0];
			$jobhash{cores} = $job->{slots};
			$jobhash{owner} = $job->{JB_owner};
			$jobhash{priority} = $job->{JB_priority};
			$jobhash{submitted} = $job->{JB_submission_time};
			$jobhash{dependency_reason} = '';
			$jobhash{children} = {};

			# Find dependency
			# If more than one (an array), assume last element is latest
			my $parents = $job->{predecessor_jobs_req};
			my $parent;
			if($parents && $parents =~ /^start_/){
				$parents = '';
			}

			if(ref($parents)){
				$parent = pop (@$parents);
			} else {
				$parent = $parents;
			}
			if($parent && length($parent) > 0){
				parse_job_dependency_hash(\%jobs, $parent, $jid, \%jobhash);
			} else {
				$jobs{$jid} = \%jobhash;
			}
		}
	}

	# Work out pipeline cwds
	# Figure out the cwd for this pipeline
	my %pipeline_wds;
	foreach my $pipelinekey (keys (%pipelines)){
		my $pipeline_wd;
		if($pipeline_single_job_ids{$pipelinekey}){
			my $cwd_command = "qstat -j $pipeline_single_job_ids{$pipelinekey} | grep cwd";
	        $pipeline_wd = `$cwd_command`;
	        $pipeline_wd =~ s/^cwd:\s+//;
	        $pipeline_wd =~ s/\n//;
			$pipeline_wds{$pipelinekey} = $pipeline_wd;
		}
	}

	# print Dumper ($data); exit;
	# print Dumper (\%jobs); exit;

	my $output = print_jobs_output(\%jobs, \%pipelines, \%pipeline_single_job_ids, \%pipeline_wds, $cols, $all_users);
	return ($output);
}


sub parse_job_dependency_hash {

	my ($hashref, $parent, $jid, $jobhash) = @_;

	foreach my $key ( keys (%{$hashref}) ){
		# If $parent is numeric, it's a job ID. If not, it's a job name
		my $jobname;
		if($parent =~ /^\d+$/){
			$jobname = $key;
		} else {
			$jobname = ${$hashref}{$key}{full_job_name};
		}
		if($jobname eq $parent){
			${$hashref}{$key}{children}{$jid} = \%$jobhash;
		} elsif (scalar(keys(%{${$hashref}{$key}{children}})) > 0){
			parse_job_dependency_hash(\%{${$hashref}{$key}{children}}, $parent, $jid, \%$jobhash);
		}
	}
}

sub print_jobs_output {

	my ($jobs, $pipelines, $pipeline_single_job_ids, $pipeline_wds, $cols, $all_users) = @_;

	# Set up counts
	my %summary_counts;
	$summary_counts{'running'} = $summary_counts{'dependency'} = $summary_counts{'pending'} = $summary_counts{'deleting'} = 0;
	$summary_counts{'pipelines'} = scalar (keys (%{$pipelines}));
	my $summary_output;
	get_job_counts(\%{$jobs}, \%summary_counts);
	if($summary_counts{'pipelines'} > 0){ 		$summary_output .= sprintf(" Cluster Flow Pipelines: %6s\n", $summary_counts{'pipelines'}); }
	if($summary_counts{'running'} > 0){ 		$summary_output .= sprintf(" Running Jobs:           %6s\n", $summary_counts{'running'}); }
	if($summary_counts{'pending'} > 0){ 		$summary_output .= sprintf(" Queued (Resources):     %6s\n", $summary_counts{'pending'}); }
	if($summary_counts{'dependency'} > 0){ 		$summary_output .= sprintf(" Queued (Dependency):    %6s\n", $summary_counts{'dependency'}); }
	if($summary_counts{'deleting'} > 0){ 		$summary_output .= sprintf(" Jobs being deleted:     %6s\n", $summary_counts{'deleting'}); }

	# Go through hash and create output
	my $output = "";
	foreach my $pipelinekey (keys (%{$pipelines})){

		my $pipeline = $pipelinekey;
		if($pipelinekey =~ /^(.+)_(\d{10})$/){
			$pipeline = $1;
		}

		# Figure out the cwd for this pipeline
		my $pipeline_wd;
		if($pipeline_wds->{$pipelinekey}){
			$pipeline_wd = $pipeline_wds->{$pipelinekey};
		}

		# See if we can find the number of jobs originally submitted
		my $sub_log = $pipeline_wd."/cf_".$pipelinekey."_submissionlog.txt";
		my $submitted_jobs;
		if(-e $sub_log){
			if(open(my $sfh, '<', $sub_log)){
				while (my $srow = <$sfh>){
					if($srow =~ /Total number of Cluster Flow jobs submitted: (\d+)/){
						$submitted_jobs = $1;
						last;
					}
				}
				close $sfh;
			}
		}
        # Get the completed jobs by subtraction
        my $completed_jobs;
        if($submitted_jobs){
            $completed_jobs = $submitted_jobs;
            $completed_jobs -= $summary_counts{$pipelinekey}{'running'} if $summary_counts{$pipelinekey}{'running'};
            $completed_jobs -= $summary_counts{$pipelinekey}{'pending'} if $summary_counts{$pipelinekey}{'pending'};
            $completed_jobs -= $summary_counts{$pipelinekey}{'dependency'} if $summary_counts{$pipelinekey}{'dependency'};
            $completed_jobs -= $summary_counts{$pipelinekey}{'deleting'} if $summary_counts{$pipelinekey}{'deleting'};
            $completed_jobs .= sprintf(" (%d%%)", ($completed_jobs/$submitted_jobs)*100);
        }

        # Make a Queued summary string
        my $queued_jobs;
        if ($summary_counts{$pipelinekey}{'pending'}){
            $queued_jobs .= $summary_counts{$pipelinekey}{'pending'}." (resources)";
            $queued_jobs .= " / " if $summary_counts{$pipelinekey}{'dependency'};
        }
        $queued_jobs .= $summary_counts{$pipelinekey}{'dependency'}." (dependencies)" if $summary_counts{$pipelinekey}{'dependency'};

		$output .= "\n".('=' x 70)."\n";
		$output .= sprintf("%-24s%-44s\n", " Cluster Flow Pipeline:", $pipeline);
		$output .= sprintf("%-24s%-44s\n", " Submitted:", CF::Helpers::parse_seconds(time - $pipelines->{$pipelinekey}{started})." ago");
		$output .= sprintf("%-24s%-44s\n", " Working Directory:", $pipeline_wd) if $pipeline_wd;
		$output .= sprintf("%-24s%-44s\n", " Cluster Flow ID:", $pipelinekey);
		$output .= sprintf("%-24s%-44s\n", " Submitted Jobs:", $submitted_jobs) if $submitted_jobs;
		$output .= sprintf("%-24s%-44s\n", " Running Jobs:", $summary_counts{$pipelinekey}{'running'}) if $summary_counts{$pipelinekey}{'running'};
		$output .= sprintf("%-24s%-44s\n", " Queued Jobs:", $queued_jobs) if $queued_jobs;
		$output .= sprintf("%-24s%-44s\n", " Deleting Jobs:", $summary_counts{$pipelinekey}{'deleting'}) if $summary_counts{$pipelinekey}{'deleting'};
		$output .= sprintf("%-24s%-44s\n", " Completed Jobs:", $completed_jobs) if $completed_jobs;
		$output .= "".('=' x 70)."\n";
		print_jobs_pipeline_output(\%{$jobs}, 0, \$output, $all_users, $cols, $pipelinekey);
	}

	# Go through jobs which don't have a pipeline
	my $notcf_output = "";
	print_jobs_pipeline_output(\%{$jobs}, 0, \$notcf_output, $all_users, $cols, 'unknown');
	if(length($notcf_output) > 0){
		$output .= "\n".('=' x 70)."\n";
		$output .= "  Not Cluster Flow Jobs  ";
		$output .= "\n".('=' x 70)."\n";
		$output .= $notcf_output;
	}

	# Go through Queing jobs which don't have a pipeline
	my $notcfpending_output = "";
	print_jobs_pipeline_output(\%{$jobs}, 0, \$notcfpending_output, $all_users, $cols, 'unknown_pending');
	if(length($notcfpending_output) > 0){
		$output .= "\n".('=' x 70)."\n";
		$output .= "  Not Cluster Flow Jobs - Queing  ";
		$output .= "\n".('=' x 70)."\n";
		$output .= $notcfpending_output;
	}

	# Print a summary of jobs - numbers are calculated at top of sub
	if($summary_output){
		$output .= "\n".('=' x 70)."\n Summary Job Counts \n".('=' x 70)."\n$summary_output";
	}

	return ("$output\n");

}

sub print_jobs_pipeline_output {

	my ($hashref, $depth, $output, $all_users, $cols, $pipeline) = @_;

	foreach my $key (keys (%{$hashref}) ){

		# Ignore this unless this is part of the pipeline we're printing
		next unless (${$hashref}{$key}{'pipelinekey'} eq $pipeline || $depth > 0);

		my $children = scalar(keys(%{${$hashref}{$key}{'children'}}));

		if($depth == 0){
			${$output} .= "\n";
		}

		${$output} .= " ".(" " x ($depth*5))."- ";

		if(${$hashref}{$key}{'state'} =~ /running/i && $cols){
			${$output} .= color 'green';
			${$output} .= " ";
		} elsif(${$hashref}{$key}{'state'} =~ /deleting/i && $cols){
			${$output} .= color 'white on_red';
			${$output} .= " ";
		} elsif ($depth == 0 && $cols) {
			${$output} .= color 'yellow';
			${$output} .= " ";
		}

		if(${$hashref}{$key}{'state'} =~ /deleting/i){
			${$output} .= "** Terminating ** ";
		}
		if(${$hashref}{$key}{'module'}){
			${$output} .= ${$hashref}{$key}{'module'};
		} else {
			${$output} .= ${$hashref}{$key}{'jobname'};
		}

		if($cols){
			${$output} .= " ";
			${$output} .= color 'reset';
		}

		# Extra info for running jobs
		# if(${$hashref}{$key}{state} eq 'running'){
		if($depth == 0){

			my @lines = split("\n", ${$output});
			my $lastline = pop(@lines);
			my $chars = length($lastline);
			my $spaces = 50 - $chars;
			${$output} .= (" " x $spaces);


			if($all_users){
				${$output} .= color 'green' if $cols;
				my $user = " {".${$hashref}{$key}{'owner'}."} ";
				${$output} .= $user;
				${$output} .= color 'reset' if $cols;
				$spaces = 12 - length($user);
				${$output} .= (" " x $spaces);
			}

			${$output} .= color 'blue' if $cols;
			my $s = ""; $s = "s" if ${$hashref}{$key}{'cores'} > 1;
			${$output} .= " [".${$hashref}{$key}{'cores'}." core$s] ";
			${$output} .= color 'reset' if $cols;

			unless(${$hashref}{$key}{'state'} =~ /running/i){
				${$output} .= color 'yellow' if $cols;
				${$output} .= " [queued, priority ".${$hashref}{$key}{'priority'}."] ";
				${$output} .= color 'reset' if $cols;
				if(length(${$hashref}{$key}{'dependency_reason'}) > 0){
					${$output} .= color 'yellow' if $cols;
					${$output} .= " (Reason: ".${$hashref}{$key}{'dependency_reason'}.")";
					${$output} .= color 'reset' if $cols;
				}
			}


			my $timestamp = "";
			my ($year, $month, $day, $hour, $minute, $second) = $timestamp =~ /^(\d{4})-(\d\d)-(\d\d)T(\d\d):(\d\d):(\d\d)/;
			if($second){
				if(${$hashref}{$key}{'started'}){
					${$output} .= color 'magenta' if $cols;
					${$output} .= "running for ";
					$timestamp = ${$hashref}{$key}{'started'};
				} else {
					${$output} .= color 'yellow' if $cols;
					${$output} .= "queued for ";
					$timestamp = ${$hashref}{$key}{'submitted'};
				}
				my $time = timelocal($second ,$minute, $hour, $day, $month-1, $year);
				my $duration = CF::Helpers::parse_seconds(time - $time, 0);
				${$output} .= $duration;
			} else {
				${$output} .= $timestamp;
			}

			${$output} .= color 'reset' if $cols;
		}
		${$output} .= "\n";

		# Now go through and print child jobs

		if($children){
			if(defined(${$hashref}{$key}{'module'}) and ${$hashref}{$key}{'module'} eq 'download'
				and defined(${$hashref}{$key}{'state'}) and ${$hashref}{$key}{'state'} ne 'running'){
				# don't increase the depth if this is a download - avoid the huge christmas trees
				print_jobs_pipeline_output(\%{${$hashref}{$key}{'children'}}, $depth, \${$output}, $all_users, $cols, $pipeline);
			} else {
				print_jobs_pipeline_output(\%{${$hashref}{$key}{'children'}}, $depth + 1, \${$output}, $all_users, $cols, $pipeline);
			}
		}
	}

	return (${$output});

}

# Built a summary of all active jobs
sub get_job_counts {
	my ($hashref, $counts) = @_;

	foreach my $key (keys (%{$hashref}) ){

        # Which pipeline is this job with?
        my $p = ${$hashref}{$key}{'pipelinekey'};

		# State for this job
		if(${$hashref}{$key}{'state'} =~ /running/i){
            ${$counts}{'running'}++;
            ${$counts}{$p}{'running'}++;
        }
		if(${$hashref}{$key}{'state'} =~ /pending/i){
			if(${$hashref}{$key}{'dependency_reason'} =~ /Dependency/i){
				${$counts}{'dependency'}++;
                ${$counts}{$p}{'dependency'}++;
			} else {
				${$counts}{'pending'}++;
                ${$counts}{$p}{'pending'}++;
			}
		}
		if(${$hashref}{$key}{'state'} =~ /deleting/i){
            ${$counts}{'deleting'}++;
            ${$counts}{$p}{'deleting'}++;
        }

		# Recursion to children if we have any
		my $children = scalar(keys(%{${$hashref}{$key}{'children'}}));
		if($children){
			get_job_counts(\%{${$hashref}{$key}{'children'}}, \%{${counts}});
		}
	}
}


# local Function to delete all jobs with a pipeline ID
sub cf_pipeline_localkill {

	my ($cf_id) = @_;
	my $kill_output;
    my $jobcount = 0;

    # Get job list
    my $osx = 0;
    my $jobs = `ps -o pid -o cmd 2>&1`;
    chomp($jobs);

    # Mac, not linux
    if($jobs =~ /cmd: keyword not found/){
        $osx = 1;
        $jobs = `ps -o pid -o command`;
        chomp($jobs);
    }

	my $output = "";
	foreach my $job (split(/\n/, $jobs)){
        # Look for top level jobs matching the id
        if($job =~ /^(\d+) bash cf_local_cf_(.+)_(\d{10})_commands.sh$/){
            my $process_id = $1;
            my $pipelinekey = "$2_$3";
            if($pipelinekey eq $cf_id){
                my $k_cmd = "kill -9 -\$(ps -o pgid= $process_id | grep -o '[0-9]*')";
                $kill_output .= `$k_cmd`;
                $jobcount++;
            }
        }
    }
	if($jobcount > 0){
		return "$jobcount cluster flow pipleline deleted.\n$kill_output\n";
	} else {
		return "Error - no jobs found for pipeline $cf_id\n";
	}

}

# SLURM Function to delete all jobs with a pipeline ID
sub cf_pipeline_scancel {

	my ($pid) = @_;
	my $jobcount = 0;
	my $scancel_output;

	# Build command
	my $curruser = `whoami`;
	my $squeue_command = 'squeue -o "%i %T %u %C %S %Q %j %E %r" -u '.$curruser;

	# Parse results
	my $squeue = `$squeue_command`;
	my @sjobs = split(/[\n\r]+/, $squeue);
	foreach my $job (@sjobs){
		# Split result
		my ($jid, $state, $owner, $cores, $started, $priority, $jobname, $dependency, $dependency_reason) = split(/ /, $job);
		next if($jid eq 'JOBID');

		my $pipelinekey;
		if($jobname =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
			$pipelinekey = "$1_$2";
		}

		if($pipelinekey && $pipelinekey eq $pid){
			my $scancel_command = "scancel $jid";
			$jobcount++;
			my $scancel = `$scancel_command`;
			$scancel_output .= $scancel;
		}
	}

	if($jobcount > 0){
		return "$jobcount jobs deleted.\n$scancel_output\n";
	} else {
		return "Error - no jobs found for pipeline $pid\n";
	}

}


# GRID Engine Function to delete all jobs with a pipeline ID
sub cf_pipeline_qdel {

	my ($pid) = @_;
	my $jobcount = 0;
	my $qdel_output;

	my $qstat_command = "qstat -pri -r -xml";
	my $qstat = `$qstat_command`;

	my $xml = new XML::Simple;
	my $data = $xml->XMLin($qstat);

	my %jobs;
	my %pipelines;
	my @jlist;

	# Running Jobs
	# Returns a hash instead of an array if only one element
	if(ref($data->{queue_info}->{job_list}) eq 'HASH'){
		@jlist = \%{$data->{queue_info}->{job_list}};
	} elsif(ref($data->{queue_info}->{job_list}) eq 'ARRAY'){
		@jlist = @{$data->{queue_info}->{job_list}};
	}
	if($data->{queue_info}->{job_list}){
		foreach my $job (@jlist){
			my $jid = $job->{JB_job_number};
			my $jobname = $job->{full_job_name};
			my $pipelinekey;

			if($jobname =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
				$pipelinekey = "$1_$2";
			}

			if($pipelinekey && $pipelinekey eq $pid){
				my $qdel_command = "qdel $jid";
				$jobcount++;
				my $qdel = `$qdel_command`;
				$qdel_output .= $qdel;
			}
		}
	}

	# Pending Jobs
	# Returns a hash instead of an array if only one element
	@jlist = ();
	if(ref($data->{job_info}->{job_list}) eq 'HASH'){
		@jlist = \%{$data->{job_info}->{job_list}};
	} elsif(ref($data->{job_info}->{job_list}) eq 'ARRAY'){
		@jlist = @{$data->{job_info}->{job_list}};
	}
	if($data->{job_info}->{job_list}){
		foreach my $job (@jlist){
			my $jid = $job->{JB_job_number};
			my $jobname = $job->{full_job_name};
			my $pipelinekey = '';

            if($jobname =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
				$pipelinekey = "$1_$2";
			}

			if($pipelinekey eq $pid){
				my $qdel_command = "qdel $jid";
				$jobcount++;
				my $qdel = `$qdel_command`;
				$qdel_output .= $qdel;
			}
		}
	}

	if($jobcount > 0){
		return "$jobcount jobs deleted:\n$qdel_output\n";
	} else {
		return "Error - no jobs found for pipeline $pid\n";
	}

}


sub cf_check_updates {

	my ($current_version) = @_;

	# Get contents of Cluster Flow current version file using LWP::Simple
	my $version_url = 'http://clusterflow.io/version.php';
	my $avail_version = get($version_url) or return ("Can't access address to check available version:\n$version_url\n\n", 0);
	my $timestamp = time();

	# Update the .cfupdates files with the available version and checked timestamp
	my $updates_file = $ENV{"HOME"}."/.clusterflow/.cfupdates";
	if(open(my $fh, '>', $updates_file)){
		# Write out the new config file contents if we can
		print($fh "$avail_version\n$timestamp\n");
		close $fh;
	}

	if(cf_compare_version_numbers($current_version, $avail_version)){
		return ("".("="x45)."\n A new version of Cluster Flow is available!\n Running v$current_version, v$avail_version available.\n".("="x45)."\n
You can download the latest version of Cluster Flow from\nhttps://github.com/ewels/clusterflow/releases/\n\n", 1);
	} else {
		return ("Your copy of Cluster Flow is up to date. Running v$current_version, v$avail_version available.\n\n", 0);
	}
}

# Function to properly split up version numbers
# Returns true if second supplied vn is newer than first
sub cf_compare_version_numbers {

	my ($vn1_string, $vn2_string) = @_;

    if(!defined($vn2_string)){
        return 0;
    }

	my @vn1_parts = split(/\.|\s+/, $vn1_string);
	my @vn2_parts = split(/\.|\s+/, $vn2_string);

	for my $i (0 .. $#vn2_parts){
		if(defined($vn1_parts[$i])){

			# Numeric checks
			if($vn1_parts[$i] =~ /^\d+$/ && $vn2_parts[$i] =~ /^\d+$/){
				if($vn2_parts[$i] > $vn1_parts[$i]){
					return 1;
				}
			} elsif($vn1_parts[$i] !~ /^\d+$/ && $vn2_parts[$i] =~ /^\d+$/){
				# 0.1.1 beats 0.1 devel
				return 1;
			}
		} else {
			return 1;
		}
	}

	return 0;

}


1; # Must return a true value
