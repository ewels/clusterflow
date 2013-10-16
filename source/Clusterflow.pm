#!/usr/bin/perl
package Clusterflow; 

use warnings;
use strict;
use FindBin qw($Bin);
use Exporter;
use Data::Dumper;

# our @ISA = qw(Exporter);
# our @EXPORT_OK = qw(load_params);



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
				# Clear out excess whitespace
				$sections[1] =~ s/^\s+//;
				$sections[1] =~ s/\s+$//;
				# Push to array
				push(@files, $sections[1]);
				$num_input_files++;
			}
		}
	}
	
	close(RUN);
	
	# If we don't have any input files, bail now
	if($num_input_files == 0){
		print "\n###  Error - no file names found from job $prev_job_id. Exiting... ###\n\n";
		exit;
	}
	
	return (\@files, $runfile, $job_id, $prev_job_id, $parameters, \%config);
}


# Function to look at supplied file names and work out whether they're paired end or not
sub is_paired_end {
	
	my @files = sort(@_);
	my @se_files;
	my @pe_files;
	
	for (my $i = 0; $i <= $#files; $i++){
		if($i < $#files){
			# Make stripped copies of the fns for comparison
			(my $fn1 = $files[$i]) =~ s/_[1-4]//;
			(my $fn2 = $files[$i+1]) =~ s/_[1-4]//;
			if($fn1 eq $fn2){
				my @pe = ($files[$i], $files[$i+1]);
				push (@pe_files, \@pe);
				$i++; # Push up $i so we ignore the next file
			} else {
				push (@se_files, $files[$i]);
			}
		} else {
			# If we have an odd number of file names and have got to the end, must be SE
			push (@se_files, $files[$i]);
		}
	}
	
	return (\@se_files, \@pe_files);
}


1; # Must return a true value