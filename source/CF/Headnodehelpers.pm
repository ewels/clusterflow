#!/usr/bin/perl
package CF::Headnodehelpers; 

use warnings;
use strict;
use FindBin qw($Bin);
use lib "$FindBin::Bin/source";
use Exporter;
use POSIX qw(strftime);
use XML::Simple;
use Time::Local;
use Term::ANSIColor;
use CF::Helpers;
use Data::Dumper;


# Function to parse qstat results and return them in a nicely formatted manner
sub parse_qstat {
	
	my ($all_users) = @_;
	
	my $qstat_command = "qstat -pri -r -xml";
	if($all_users){
		$qstat_command .= ' -u "*"';
	}
	
	my $qstat = `$qstat_command`;
	
	my $xml = new XML::Simple;
	my $data = $xml->XMLin($qstat);
	
	# print Dumper $data; exit;
	
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
			my $pipeline = 'unknown';
			my $jobname = $job->{full_job_name};
			
			if($jobname =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
				$pipeline = $1;
				$pipelines{$pipeline}{started} = $2;
				$jobs{$jid}{module} = $3;
			}
			
			$jobs{$jid}{pipeline} = $pipeline;
			$jobs{$jid}{jobname} = $jobname;
			$jobs{$jid}{state} =  $job->{state}->[0];
			if($job->{state}->[1] eq 'dr'){
				$jobs{$jid}{state} = 'deleting';
			}
			$jobs{$jid}{cores} = $job->{slots};
			$jobs{$jid}{mem} = $job->{hard_request}->{content};
			$jobs{$jid}{owner} = $job->{JB_owner};
			$jobs{$jid}{priority} = $job->{JB_priority};
			$jobs{$jid}{started} = $job->{JAT_start_time};
			$jobs{$jid}{children} = {};
			
		}
	}
	
	# Pending Jobs
	foreach my $job (@{$data->{job_info}->{job_list}}){
		my $jid = $job->{JB_job_number};
		my %jobhash;
		my $pipeline = 'unknown_pending';
		my $jobname = $job->{full_job_name};
		
		if($jobname =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
			$pipeline = $1;
			$pipelines{$pipeline}{started} = $2;
			$jobhash{module} = $3;
		}
		
		$jobhash{pipeline} = $pipeline;
		$jobhash{jobname} = $jobname;
		$jobhash{state} = $job->{state}->[0];
		$jobhash{cores} = $job->{slots};
		$jobhash{owner} = $job->{JB_owner};
		$jobhash{priority} = $job->{JB_priority};
		$jobhash{submitted} = $job->{JB_submission_time};
		$jobhash{children} = {};
		
		# Find dependency
		# If more than one (an array), assume last element is latest
		my $parents = $job->{predecessor_jobs_req};
		my $parent;
		if(ref($parents)){
			$parent = pop (@$parents);
		} else {
			$parent = $parents;
		}
		if($parent && length($parent) > 0){
			parse_qstat_search_hash(\%jobs, $parent, $jid, \%jobhash);
		} else {
			$jobs{$jid} = \%jobhash;
		}
	}
	
	# print Dumper ($data); exit;
	# print Dumper (\%jobs); exit;
	
	# Go through hash and create output
	my $output = "";
	foreach my $pipeline (keys (%pipelines)){
		$output .= "\n".('=' x 50)."\n";
		$output .= "  Cluster Flow Pipeline $pipeline\n";
		$output .= "  Submitted ".CF::Helpers::parse_seconds(time - $pipelines{$pipeline}{started})." ago";
		$output .= "\n".('=' x 50)."\n";
		parse_qstat_print_hash(\%jobs, 0, \$output, $all_users, $pipeline);
	}
	
	# Go through jobs which don't have a pipeline
	my $notcf_output = "";
	parse_qstat_print_hash(\%jobs, 0, \$notcf_output, $all_users, 'unknown');
	if(length($notcf_output) > 0){
		$output .= "\n".('=' x 50)."\n";
		$output .= "  Not Cluster Flow Jobs  ";
		$output .= "\n".('=' x 50)."\n";
		$output .= $notcf_output;
	}
	
	my $notcfpending_output = "";
	parse_qstat_print_hash(\%jobs, 0, \$notcfpending_output, $all_users, 'unknown_pending');
	if(length($notcfpending_output) > 0){
		$output .= "\n".('=' x 50)."\n";
		$output .= "  Not Cluster Flow Jobs - Queing  ";
		$output .= "\n".('=' x 50)."\n";
		$output .= $notcfpending_output;
	}
	
	return ("$output\n");	
	
}

sub parse_qstat_search_hash {

	my ($hashref, $parent, $jid, $jobhash) = @_;
	
	foreach my $key (keys (%{$hashref}) ){
		my $jobname = ${$hashref}{$key}{jobname};
		if($jobname eq $parent){
			${$hashref}{$key}{children}{$jid} = \%$jobhash;
		} elsif (scalar(keys(%{${$hashref}{$key}{children}})) > 0){
			foreach my $child (keys %{${$hashref}{$key}{children}}){
				parse_qstat_search_hash(\%{${$hashref}{$key}{children}}, $parent, $jid, $jobhash);
			}
		}
	}
}

sub parse_qstat_print_hash {

	my ($hashref, $depth, $output, $all_users, $pipeline) = @_;

	foreach my $key (keys (%{$hashref}) ){
	
		# Ignore this unless this is part of the pipeline we're printing
		next unless (${$hashref}{$key}{pipeline} eq $pipeline || $depth > 0);
		
		my $children = scalar(keys(%{${$hashref}{$key}{children}}));
		
		if($depth == 0){
			${$output} .= "\n";
		}
		
		${$output} .= " ".(" " x ($depth*5))."- ";
		
		if(${$hashref}{$key}{state} eq 'running'){
			${$output} .= color 'red on_white';
			${$output} .= " ";
		} elsif(${$hashref}{$key}{state} eq 'deleting'){
			${$output} .= color 'white on_red';
			${$output} .= " ";
		} elsif ($depth == 0) {
			${$output} .= color 'yellow on_white';
			${$output} .= " ";
		}
		if(${$hashref}{$key}{module}){
			${$output} .= ${$hashref}{$key}{module};
		} else {
			${$output} .= ${$hashref}{$key}{jobname};
		}
		${$output} .= " ";
		${$output} .= color 'reset';
		
		# Extra info for running jobs
		# if(${$hashref}{$key}{state} eq 'running'){
		if($depth == 0){
			
			my @lines = split("\n", ${$output});
			my $lastline = pop(@lines);
			my $chars = length($lastline);
			my $spaces = 50 - $chars;
			${$output} .= (" " x $spaces);
			
			
			if($all_users){
				${$output} .= color 'green';
				my $user = " {".${$hashref}{$key}{owner}."} ";
				${$output} .= $user;
				${$output} .= color 'reset';
				$spaces = 12 - length($user);
				${$output} .= (" " x $spaces);
			}
			
			${$output} .= color 'blue';
			${$output} .= " [".${$hashref}{$key}{cores}." cores] ";
			${$output} .= color 'reset';
			
			if(${$hashref}{$key}{state} ne 'running'){
				${$output} .= color 'yellow';
				${$output} .= " [priority ".${$hashref}{$key}{priority}."] ";
				${$output} .= color 'reset';
			}
			
			
			my $timestamp = "";
			if(${$hashref}{$key}{started}){
				${$output} .= color 'magenta';
				${$output} .= "running for ";
				$timestamp = ${$hashref}{$key}{started};
			} else {
				${$output} .= color 'yellow';
				${$output} .= "queued for ";
				$timestamp = ${$hashref}{$key}{submitted};
			}
			my ($year, $month, $day, $hour, $minute, $second) = $timestamp =~ /^(\d{4})-(\d\d)-(\d\d)T(\d\d):(\d\d):(\d\d)/;
			my $time = timelocal($second ,$minute, $hour, $day, $month-1, $year);
			my $duration = CF::Helpers::parse_seconds(time - $time, 0);
			${$output} .= $duration;
			
			${$output} .= color 'reset';
		}
		${$output} .= "\n";
		
		# Now go through and print child jobs
		if($children){
			foreach my $child (keys %{${$hashref}{$key}{children}}){
				parse_qstat_print_hash(\%{${$hashref}{$key}{children}}, $depth + 1, \${$output}, $all_users, $pipeline);
			}
		}
	}
	
	return (${$output});

}




1; # Must return a true value