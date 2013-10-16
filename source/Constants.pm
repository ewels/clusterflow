#!/usr/bin/perl
package Constants;

use warnings;
use strict;
use FindBin;



# This is the global version for Cluster Flow; it is updated
# manually in new releases.  You don't need to (and
# shouldn't) edit this yourself

our $CF_VERSION = "0.1 devel";

# our $homedir = File::HomeDir->my_home;
our $homedir = $ENV{"HOME"};

# Old config hash. Delete soon.
our %config;

# Empty config vars
our $EMAIL;
our @NOTIFICATIONS;

# Hard coded defaults
our $PRIORITY = -500;
our $TOTAL_CORES = 64;
our $DEFAULT_MEM = '4G';








parse_conf_file ();

sub parse_conf_file {
	
	
	# Read global config variables in. Do in order so that local prefs overwrite.
	
	my @config_files = ("$FindBin::Bin/clusterflow.config", "$homedir/clusterflow/clusterflow.config", './clusterflow.config');
	foreach my $config_file (@config_files){
		if(-e $config_file){
			open (CONFIG, $config_file) or die "Can't read $config_file: $!";
			my $comment_block = 0;
			while (<CONFIG>) {
				chomp;
				s/\n//;
				s/\r//;
				if($_ =~ /^\/\*/){
					$comment_block = 1;
					next;
				}
				if($_ =~ /^\*\//){
					$comment_block = 0;
					next;
				}
				if($_ =~ /^\@/ && !$comment_block){
					my @sections = split(/:/, $_, 2);
					$config{substr($sections[0], 1)} = $sections[1];
					my $name = substr($sections[0], 1);
					my $val = $sections[1];
					
					if($name eq 'email'){
						$EMAIL = $val;
					} elsif($name eq 'notification'){
						push @NOTIFICATIONS, $val;
					} elsif($name eq 'priority'){
						$PRIORITY = $val;
					} elsif($name eq 'total_cores'){
						$TOTAL_CORES = $val;
					} elsif($name eq 'default_mem'){
						$DEFAULT_MEM = $val;
					}
				}
			}
			close(CONFIG);
		}
	}

}


1;