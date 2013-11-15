#!/usr/bin/perl
use warnings;
use strict;
use FindBin qw($Bin);
use Exporter;

package CF::Constants;

our $CF_VERSION = "0.1 devel";

# our $homedir = File::HomeDir->my_home;
our $homedir = $ENV{"HOME"};

# Old config hash. Delete soon.
our %config;

# Empty config vars
our $EMAIL;
our @NOTIFICATIONS;
our %GENOME_PATHS;
our %BOWTIE_PATHS;
our %GTF_PATHS;

# Hard coded defaults
our $SPLIT_FILES = 1;
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
				
				if($_ =~ /^\/\*.*\*\/$/){	# one line comments
					$comment_block = 0;
					next;
				} elsif($_ =~ /^\/\*/){		# multiline comments start
					$comment_block = 1;
					next;
				} elsif($_ =~ /^\*\//){		# multiline comments start
					$comment_block = 0;
					next;
				}
				if($_ =~ /^\@/ && !$comment_block){
					my @sections = split(/\t/, $_);
					$config{substr($sections[0], 1)} = $sections[1];
					my $name = substr($sections[0], 1);
					my $val = $sections[1];
					
					if($name eq 'email'){
						$EMAIL = $val;
					} elsif($name eq 'notification'){
						push @NOTIFICATIONS, $val;
					} elsif($name eq 'genome_path'){
						$GENOME_PATHS{$val} = $sections[2];
					} elsif($name eq 'bowtie_path'){
						$BOWTIE_PATHS{$val} = $sections[2];
					} elsif($name eq 'gtf_path'){
						$GTF_PATHS{$val} = $sections[2];
					} elsif($name eq 'split_files'){
						$SPLIT_FILES = $val;
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
	
	# Remove duplicate Notifications
	my @unique_notifications;
	my %seen_notification;
	foreach my $value (@NOTIFICATIONS) {
		if (!$seen_notification{$value}) {
			push @unique_notifications, $value;
			$seen_notification{$value} = 1;
		}
	}
	@NOTIFICATIONS = @unique_notifications;
}


sub runfile_constants {

	my $output;
	
	$output = <<"EOT";
\@email	$EMAIL
\@split_files	$SPLIT_FILES
\@priority	$PRIORITY
\@total_cores	$TOTAL_CORES
\@default_mem	$DEFAULT_MEM
EOT
	
	foreach my $not (@NOTIFICATIONS){
		$output .= "\@notification\t$not\n";
	}
		
	return ($output);
	
}

# Prints help for a specific module or pipeline
sub clusterflow_pipeline_help {
	
	my ($pipeline) = @_;
	
	my $help = "";
	
	my @pipelines = ('./$pipeline.config', "$homedir/clusterflow/pipelines/$pipeline.config", "$FindBin::Bin/pipelines/$pipeline.config");
	my @modules = ("$homedir/clusterflow/modules/$pipeline", "$FindBin::Bin/modules/$pipeline");
	foreach my $pipeline (@pipelines){
		if(-e $pipeline){
			open (PIPELINE, $pipeline) or die "Can't read $pipeline: $!";
			my $comment_block = 0;
			while (<PIPELINE>) {
				chomp;
				s/\n//;
				s/\r//;
				if($_ =~ /^\/\*/){		# multiline comments start
					$comment_block = 1;
				} elsif($_ =~ /^\*\//){		# multiline comments start
					$comment_block = 0;
				} elsif($comment_block){
					$help .= $_."\n";
				}
			}
			close(PIPELINE);
		}
		if($help){
			return ($help);
		}
	}
	
	foreach my $module (@modules){
		if(-e $module){
			$help = `$module --help`;
			return ($help);
		}
	}
	
	if($help eq ""){
		$help = "\nSorry, no help found for this pipeline.\n\n";
	}
	
	return ($help);
	
}

# Prints main cluster flow help
sub clusterflow_help {

	my $help;
	
	$help = <<"EOT";

Cluster Flow Help
=================
Running Cluster Flow version $CF_VERSION

SYNTAX
	cf [flags] pipeline_name file_1 file_2..

EXAMPLE
	cf --genome NCBIM37 sra_bismark *.sra

SPECIFIC PIPELINE / MODULE HELP
	To see specific help about a pipeline or module, use
	cf --help followed by a pipeline or module name.

INTRODUCTION
	Cluster Flow is simple package to run pipelines in a cluster environment.
	
	Cluster Flow will set off multiple queued jobs on the cluster with queue 
	dependencies as defined in the pipeline.

AVAILABLE FLAGS
	For a full description of the avilable flags and how to use them, see
	the Cluster Flow documentation.

	--genome <ID>
		ID of a genome referred to in clusterflow.config
		This genome ID is used to specify genome paths, bowtie
		index basenames and GTF file paths.
		
	--genome_path <path>
		Path to a genome to be used for alignment. Overrides --genome
		
	--bowtie_path <path>
		Path to a bowtie index basename to be used for alignment.
		Overrides any bowtie path set with --genome
		
	--gtf_path <path>
		Path to a GTF file to be used for alignment (eg. Tophat).
		Overrides any GTF path set with --genome
		
	--paired
		Force paired-end mode
		
	--single
		Force single-end mode
		
	--file_list
		Text file containing input files or download URLs
		
	--module
		Run a single module instead of a pipeline
		
	--split-files <num>
		Create one run per <num> files
		
	--email <email>
		Set the e-mail address for notifications
		
	--priority <num>
		Set the queue priority for cluster jobs
		
	--cores <num>
		Set the maximum number of cores to use for all runs
		
	--mem <string>
		Set the maximum memory to use for all runs
		
	--notifications <cresa>
		Specify desired notifications
		c = pipeline complete, r = run complete, e = qsub job ends
		s = qsub job suspended, a = qsub job aborted
		
	--list_pipelines
		Print available pipelines
		
	--list_modules
		Print available modules
		
	--list_genomes
		Print available genomes
		
	--dryrun
		Prints jobs to terminal instead of submitting them to the cluster
	
	--make_config
		Interactive prompt to generate a personalised CF config file
	
	--version
		Print version of Cluster Flow installed
		
	--help
		Print this help message


AUTHOR
	Written by Phil Ewels, Babraham Institute.
	
SEE ALSO
	There is a full Cluster Flow manual available.

EOT
	
	return ($help);

}




# Function to run interactive shell prompt to generate a config file for first run
sub clusterflow_make_config{
	
	my $fn = $homedir."/clusterflow/clusterflow.config";
	
	print "\n\nCluster Flow Config Generator\n==============================\nRunning Cluster Flow version $CF_VERSION\n";
	print "This mode will generate a personalised Cluster Flow config file for you \nin your home directory: ~/clusterflow/clusterflow.config\n\n";
	
	if(-e $fn){
		print "### WARNING ###\n$fn already exists!\nThis script will overwrite that file. Do you want to continue? (y/n)\n\n";
		while (my $continue = <STDIN>){
			chomp ($continue);
			if ($continue =~ /^n(o)?/i){
				print "\nProbably wise.. See the manual for more information about how the config file works.\n\n";
				exit;
			} elsif($continue =~ /^y(es)?/i){
				print "\nOk, no problem.. I'll wipe it when we get to the end.\n\n";
				last;
			} else {
				print "\nSorry, I didn't understand that.\nCould you try again please? (y/n)\n\n";
			}
		}
	}
	
	my $config = "/*
Clusterflow Config
-------------------
Default static variables for clusterflow.
Syntax - \@key:value
These will overwrite any with the same name in the centralised config file
-------------------
*/\n\n";
	
	my $email;
	print "Right, let's get started! First off - what is your e-mail address?\nThis will be used for sending notifications.\n\n";
	while ($email = <STDIN>){
		chomp ($email);
		if($email =~ /^\w[\w\.\-]*\w\@\w[\w\.\-]*\w(\.\w{2,4})$/){
			print "\nGreat! That looks good..\n\n";
			last;
		} else {
			print "\nHmm, that e-mail address looked a little odd.\nAre you sure you typed it in correctly? (y/n)\n\n";
			my $invalidemail = <STDIN>;
			chomp($invalidemail);
			if($invalidemail =~ /^y(es)?/i){
				print "\nFair enough!\n\n";
				last;
			} else {
				print "\nNo problem.. Please try it again..\n\n";
			}
		}
	}
	$config .= "\@email	$email\n";
	
	my ($notify_complete, $notify_run, $notify_success, $notify_error, $notify_abort);
	
	print "Would you like to receive a notification when a pipeline is completed?
The e-mail tells you the pipeline that has finished, the working directory
for that pipeline, a list of Cluster Flow highlight notifications (typically
whether each step in the pipeline ran successfully for each file) and then the
log file output for each file. These notifications are recommended (y/n)\n\n";
	while ($notify_complete = <STDIN>){
		chomp ($notify_complete);
		if($notify_complete =~ /^y(es)?/i){
			print "\nGreat!\n\n";
			$config .= "\@notification	complete\n";
			last;
		} elsif ($notify_complete =~ /^n(o)?/i){
			print "\nOk, fair enough..\n\n";
			last;
		} else {
			print "\nSorry, I didn't understand that.\nCould you try again please? (y/n)\n\n";
		}
	}
	
	print "Would you like to receive a notification when each run is completed?
Usually a run is the processing of one input file. The e-mail tells you the
name of the run that has finished, its pipeline and the working directory.
It includes a list of Cluster Flow highlight notifications (typically
whether each step in the pipeline ran successfully) and then the log file output.
These notifications are recommended for those who like to keep a close eye on
their processing (y/n)\n\n";
	while ($notify_run = <STDIN>){
		chomp ($notify_run);
		if($notify_run =~ /^y(es)?/i){
			print "\nGreat!\n\n";
			$config .= "\@notification	run\n";
			last;
		} elsif ($notify_run =~ /^n(o)?/i){
			print "\nOk, sounds good..\n\n";
			last;
		} else {
			print "\nSorry, I didn't understand that.\nCould you try again please? (y/n)\n\n";
		}
	}
	
	print "Would you like to receive a notification when step of each run ends?
This will be a GRID Engine notice for every qsub job. These notifications
are not recommended as a typicaly Cluster Flow run can flood your inbox with hundreds
of such e-mails. Would you like to receive them? (y/n)\n\n";
	while ($notify_success = <STDIN>){
		chomp ($notify_success);
		if($notify_success =~ /^y(es)?/i){
			print "\nFair enough, you were warned!\n\n";
			$config .= "\@notification	end\n";
			last;
		} elsif ($notify_success =~ /^n(o)?/i){
			print "\nProbably sensible..\n\n";
			last;
		} else {
			print "\nSorry, I didn't understand that.\nCould you try again please? (y/n)\n\n";
		}
	}
	
	print "Would you like to receive a notification when a GRID Engine
job is suspended? You're unlikely to get many if any, so they're recommended.
Would you like to receive these notifications? (y/n)\n\n";
	while ($notify_error = <STDIN>){
		chomp ($notify_error);
		if($notify_error =~ /^y(es)?/i){
			print "\nSounds good!\n\n";
			$config .= "\@notification	suspend\n";
			last;
		} elsif ($notify_error =~ /^n(o)?/i){
			print "\nFair enough..\n\n";
			last;
		} else {
			print "\nSorry, I didn't understand that.\nCould you try again please? (y/n)\n\n";
		}
	}
	
	print "Ok, last one. Would you like to receive a notification 
when a GRID Engine job exits in an abort state? This typically only
happens when you or an administrator kills your cluster jobs using
qdel. You're unlikely to get many of these, so they're recommended.
Would you like to receive these notifications? (y/n)\n\n";
	while ($notify_abort = <STDIN>){
		chomp ($notify_abort);
		if($notify_abort =~ /^y(es)?/i){
			print "\nSounds good!\n\n";
			$config .= "\@notification	abort\n";
			last;
		} elsif ($notify_abort =~ /^n(o)?/i){
			print "\nFair enough..\n\n";
			last;
		} else {
			print "\nSorry, I didn't understand that.\nCould you try again please? (y/n)\n\n";
		}
	}
	$config .= "\n\n\n";
	
	
	print "\n\nGreat, that's it! The following config file will be created:\n\n$config\n";
	
	print "\nRemember that you can add further settings to your
personalised config file - see the Cluster Flow manual
for further information.\n\n\n";
	
	unless(-e $homedir."/clusterflow/" && -d $homedir."/clusterflow/"){
		mkdir ($homedir."/clusterflow/") or die "Can't create clusterflow directory: $!";
	}
	open (OUT, '>', $fn) or die "Can't write to $fn: $!";
	print OUT $config;
	close OUT;
	

}


1;