#!/usr/bin/perl
use warnings;
use strict;
use FindBin qw($Bin);
use Exporter;

# our @ISA = qw(Exporter);
# our @EXPORT_OK = qw(load_params);

package Clusterflow;

sub load_runfile_params {
	my ($runfile, $job_id, $prev_job_id, $parameters) = @_;
	unless (defined $prev_job_id && length($prev_job_id) > 0) {
		die "Previous job ID not specified\n";
	}
	unless ($parameters) {
		$parameters = '';
	}
	
	my $num_input_files = 0;
	
	warn "\nRun File:\t\t$runfile\nJob ID:\t\t\t$job_id\nPrevious Job ID:\t$prev_job_id\nParameters:\t\t$parameters\n\n";

	open (RUN,$runfile) or die "Can't read $runfile: $!";

	my @files;
	my %config;
	my $comment_block = 0;
	while(<RUN>){
	
		# clean up line
		chomp;
		s/\n//;
		s/\r//;
		
		# Ignore comment blocks
		if($_ =~ /^\/\*/){
			$comment_block = 1;
			next;
		}
		if($_ =~ /^\*\//){
			$comment_block = 0;
			next;
		}
		
		# Get config variables
		if($_ =~ /^\@/ && !$comment_block){
			my @sections = split(/:/, $_, 2);
			$config{substr($sections[0], 1)} = $sections[1];
		}
		
		# Get files
		if($_ =~ /^[^@#]/ && !$comment_block){
			my @sections = split(/\t/, $_, 2);
			if($sections[0] eq $prev_job_id){
				push(@files, $sections[1]);
				$num_input_files++;
			}
		}
	}
	
	close(RUN);
	
	# If we don't have any input files, bail now
	if($num_input_files == 0){
		print "Error - no input files found.\n";
		exit;
	}
	
	return (\@files, $runfile, $job_id, $prev_job_id, $parameters, \%config);
}




1; # Must return a true value