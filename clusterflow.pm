#!perl
use warnings;
use strict;
use Exporter;

package clusterflow;

sub load_params {
	my ($runfile, $job_id, $prev_id, $parameters) = @ARGV;
	unless (defined $prev_id && length($prev_id) > 0) {
		die "Previous job ID not specified\n";
	}

	open (RUN,$runfile) or die "Can't read $runfile: $!";

	$my @files;
	while(<>){
		chomp;
		my @sections = split(/\t/, $_, 2);
		if($sections[0] eq $prev_id){
			push(@files, $sections[1]);
		}
	}
	
	close(RUN);
	
	open (RUN,'>>',$runfile) or die "Can't write to $runfile: $!";
}

sub close_module {
	close (RUN);
}


1; # Must return a true value