#!perl
use warnings;
use strict;
use Exporter;

# our @ISA = qw(Exporter);
# our @EXPORT_OK = qw(load_params);

package Clusterflow;

sub load_params {
	my ($runfile, $job_id, $prev_job_id, $parameters) = @_;
	unless (defined $prev_job_id && length($prev_job_id) > 0) {
		die "Previous job ID not specified\n";
	}
	
	warn "\nRun File:\t\t$runfile\nJob ID:\t\t\t$job_id\nPrevious Job ID:\t$prev_job_id\nParameters:\t\t$parameters\n\n";

	open (RUN,$runfile) or die "Can't read $runfile: $!";

	my @files;
	my $comment_block = 0;
	while(<RUN>){
		chomp;
		if($_ =~ /^\/\*/){
			$comment_block = 1;
			next;
		}
		if($_ =~ /^\*\//){
			$comment_block = 0;
			next;
		}
		if($_ =~ /^[^@#]/ && !$comment_block){
			my @sections = split(/\t/, $_, 2);
			if($sections[0] eq $prev_job_id){
				push(@files, $sections[1]);
			}
		}
	}
	
	close(RUN);
	
	return (\@files, $runfile, $job_id, $prev_job_id, $parameters);
}




1; # Must return a true value