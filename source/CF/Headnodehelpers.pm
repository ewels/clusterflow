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
use LWP::Simple;

# Function to parse qstat results and return them in a nicely formatted manner
sub parse_qstat {
	
	my ($all_users, $cols) = @_;
	
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
			my $pipelinekey = 'unknown';
			my $jobname = $job->{full_job_name};
			
			if($jobname =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
				$pipeline = $1;
				$pipelinekey = "$1_$2";
				$pipelines{$pipelinekey}{started} = $2;
				$jobs{$jid}{module} = $3;
			}
			
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
			
			$jobhash{pipeline} = $pipeline;
			$jobhash{pipelinekey} = $pipelinekey;
			$jobhash{jobname} = $jobname;
			$jobhash{full_job_name} = $job->{full_job_name};
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
	}
	
	# print Dumper ($data); exit;
	# print Dumper (\%jobs); exit;
	
	# Go through hash and create output
	my $output = "";
	foreach my $pipelinekey (keys (%pipelines)){
		my $pipeline = $pipelinekey;
		if($pipelinekey =~ /^(.+)_(\d{10})$/){
			$pipeline = $1;
		}
		$output .= "\n".('=' x 50)."\n";
		$output .= "  Cluster Flow Pipeline $pipeline\n";
		$output .= "  Submitted ".CF::Helpers::parse_seconds(time - $pipelines{$pipelinekey}{started})." ago";
		$output .= "\n".('=' x 50)."\n";
		parse_qstat_print_hash(\%jobs, 0, \$output, $all_users, $cols, $pipelinekey);
	}
	
	# Go through jobs which don't have a pipeline
	my $notcf_output = "";
	parse_qstat_print_hash(\%jobs, 0, \$notcf_output, $all_users, $cols, 'unknown');
	if(length($notcf_output) > 0){
		$output .= "\n".('=' x 50)."\n";
		$output .= "  Not Cluster Flow Jobs  ";
		$output .= "\n".('=' x 50)."\n";
		$output .= $notcf_output;
	}
	
	# Go through Queing jobs which don't have a pipeline
	my $notcfpending_output = "";
	parse_qstat_print_hash(\%jobs, 0, \$notcfpending_output, $all_users, $cols, 'unknown_pending');
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
	
	foreach my $key ( keys (%{$hashref}) ){
		my $jobname = ${$hashref}{$key}{full_job_name};
		if($jobname eq $parent){
			${$hashref}{$key}{children}{$jid} = \%$jobhash;
		} elsif (scalar(keys(%{${$hashref}{$key}{children}})) > 0){
			parse_qstat_search_hash(\%{${$hashref}{$key}{children}}, $parent, $jid, \%$jobhash);
		}
	}
}

sub parse_qstat_print_hash {

	my ($hashref, $depth, $output, $all_users, $cols, $pipeline) = @_;

	foreach my $key (keys (%{$hashref}) ){
	
		# Ignore this unless this is part of the pipeline we're printing
		next unless (${$hashref}{$key}{pipelinekey} eq $pipeline || $depth > 0);
		
		my $children = scalar(keys(%{${$hashref}{$key}{children}}));
		
		if($depth == 0){
			${$output} .= "\n";
		}
		
		${$output} .= " ".(" " x ($depth*5))."- ";
		
		if(${$hashref}{$key}{state} eq 'running' && $cols){
			${$output} .= color 'red on_white';
			${$output} .= " ";
		} elsif(${$hashref}{$key}{state} eq 'deleting' && $cols){
			${$output} .= color 'white on_red';
			${$output} .= " ";
		} elsif ($depth == 0 && $cols) {
			${$output} .= color 'yellow on_white';
			${$output} .= " ";
		}
		
		if(${$hashref}{$key}{state} eq 'deleting'){
			${$output} .= "** Terminating ** ";
		}
		if(${$hashref}{$key}{module}){
			${$output} .= ${$hashref}{$key}{module};
		} else {
			${$output} .= ${$hashref}{$key}{jobname};
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
				my $user = " {".${$hashref}{$key}{owner}."} ";
				${$output} .= $user;
				${$output} .= color 'reset' if $cols;
				$spaces = 12 - length($user);
				${$output} .= (" " x $spaces);
			}
			
			${$output} .= color 'blue' if $cols;
			${$output} .= " [".${$hashref}{$key}{cores}." cores] ";
			${$output} .= color 'reset' if $cols;
			
			if(${$hashref}{$key}{state} ne 'running'){
				${$output} .= color 'yellow' if $cols;
				${$output} .= " [priority ".${$hashref}{$key}{priority}."] ";
				${$output} .= color 'reset' if $cols;
			}
			
			
			my $timestamp = "";
			if(${$hashref}{$key}{started}){
				${$output} .= color 'magenta' if $cols;
				${$output} .= "running for ";
				$timestamp = ${$hashref}{$key}{started};
			} else {
				${$output} .= color 'yellow' if $cols;
				${$output} .= "queued for ";
				$timestamp = ${$hashref}{$key}{submitted};
			}
			my ($year, $month, $day, $hour, $minute, $second) = $timestamp =~ /^(\d{4})-(\d\d)-(\d\d)T(\d\d):(\d\d):(\d\d)/;
			my $time = timelocal($second ,$minute, $hour, $day, $month-1, $year);
			my $duration = CF::Helpers::parse_seconds(time - $time, 0);
			${$output} .= $duration;
			
			${$output} .= color 'reset' if $cols;
		}
		${$output} .= "\n";
		
		# Now go through and print child jobs
		if($children){
			parse_qstat_print_hash(\%{${$hashref}{$key}{children}}, $depth + 1, \${$output}, $all_users, $cols, $pipeline);
		}
	}
	
	return (${$output});

}


sub cf_check_updates {

	my ($current_version) = @_;
	$current_version =~ s/[^\d.]//g;

	
	# Get contents of Cluster Flow current version file using LWP::Simple
	# my $version_url = 'http://www.bioinformatics.babraham.ac.uk/projects/cluster_flow/version.txt';
	my $version_url = 'http://bilin1/projects/cluster_flow/version.txt';
	my $avail_version = get($version_url) or die "Unable to fetch version: $version_url\n";
	
	$avail_version =~ s/[^\d.]//g;
	
	# Update the config files with the available version
	my @config_files = ("$FindBin::Bin/clusterflow.config", $ENV{"HOME"}."/clusterflow/clusterflow.config", './clusterflow.config');
	foreach my $config_file (@config_files){
		if(-e $config_file){
			
			# Read in the config file contents
			my $config_file_contents;
			{
				open (my $fh, $config_file) or die "Can't read $config_file: $!";
				local $/;
				$config_file_contents = <$fh>;
				close $fh;
			}
			
			# Swap existing variables with new ones
			my $timestamp = time();
			$config_file_contents =~ s/^\@available_version(.*)\n*/\@available_version\t$avail_version\n/img;
			$config_file_contents =~ s/^\@updates_last_checked(.*)\n*/\@updates_last_checked\t$timestamp\n/img;
			
			# If we don't have the variables, add them
			if($config_file_contents !~ /\@available_version/){
				$config_file_contents .= "\n\@available_version\t$avail_version\n";
			}
			if($config_file_contents !~ /\@updates_last_checked/){
				$config_file_contents .= "\n\@updates_last_checked\t$timestamp\n";
			}
			
			# Write out the new config file contents
			{
				open(my $fh, '>', $config_file) or die "Can't create $config_file: $!\n";
				print($fh $config_file_contents);
				close $fh;
			}
		}
	}
	
	if($avail_version > $current_version){
		return "".("="x45)."\n A new version of Cluster Flow is available!\n Running v$current_version, v$avail_version available.\n".("="x45)."\n
You can download the latest version of Cluster Flow from\nhttp://www.bioinformatics.babraham.ac.uk/projects/cluster_flow/\n\n";
	} else {
		return "Your copy of Cluster Flow is up to date. Running v$current_version, v$avail_version available.\n\n";
	}
}



1; # Must return a true value