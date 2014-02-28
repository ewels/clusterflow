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
our $CHECK_UPDATES;
our @NOTIFICATIONS;
our $SPLIT_FILES = 1;
our $PRIORITY = -500;
our $TOTAL_CORES = 64;
our $TOTAL_MEM = '4G';
our $CLUSTER_ENVIRONMENT = 'GRIDEngine';
our $CF_MODULES = 1;

# Empty genome path vars
our %GENOME_PATH_CONFIGS;
our %BOWTIE_PATH_CONFIGS;
our %GTF_PATH_CONFIGS;
our %GENOME_PATHS;
our %BOWTIE_PATHS;
our %GTF_PATHS;
our %GENOME_SPECIES;
our %BOWTIE_SPECIES;
our %GTF_SPECIES;
our %GENOME_ASSEMBLIES;
our %BOWTIE_ASSEMBLIES;
our %GTF_ASSEMBLIES;

# Update checking variables
our $AVAILABLE_VERSION;
our $UPDATES_LAST_CHECKED = 0;

parse_conf_file ();
parse_genomes_file();
parse_updates_file();

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
					my @sections = split(/\s+/, $_);
					$config{substr($sections[0], 1)} = $sections[1];
					my $name = substr($sections[0], 1);
					my $val = $sections[1];
					
					if($name eq 'email'){
						$EMAIL = $val;
					} elsif($name eq 'check_updates'){
						$CHECK_UPDATES = $val;
					} elsif($name eq 'available_version'){
						$AVAILABLE_VERSION = $val;
					} elsif($name eq 'updates_last_checked'){
						$UPDATES_LAST_CHECKED = $val;
					} elsif($name eq 'notification'){
						push @NOTIFICATIONS, $val;
					} elsif($name eq 'split_files'){
						$SPLIT_FILES = $val;
					} elsif($name eq 'priority'){
						$PRIORITY = $val;
					} elsif($name eq 'total_cores'){
						$TOTAL_CORES = $val;
					} elsif($name eq 'total_mem'){
						$TOTAL_MEM = $val;
					} elsif($name eq 'cluster_environment'){
						$CLUSTER_ENVIRONMENT = $val;
					} elsif($name eq 'ignore_modules'){
						$CF_MODULES = 0;
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

sub parse_genomes_file {
	
	# Read genomes config variables in. Do in order so that local prefs overwrite.
	
	my @genome_files = ("$FindBin::Bin/genomes.config", "$homedir/clusterflow/genomes.config", './genomes.config');
	foreach my $genome_file (@genome_files){
		if(-e $genome_file){
			open (GCONFIG, $genome_file) or die "Can't read $genome_file: $!";
			my $comment_block = 0;
			while (<GCONFIG>) {
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
					my @sections = split(/\s+/, $_);
					$config{substr($sections[0], 1)} = $sections[1];
					my $path_type = substr($sections[0], 1);
					my $key = $sections[1];
					my $path = $sections[2];
					my $species = $sections[3];
					my $assembly = $sections[4];
					if($path_type eq 'genome_path'){
						$GENOME_PATH_CONFIGS{$key} = $genome_file;
						$GENOME_PATHS{$key} = $path;
						$GENOME_SPECIES{$key} = $species;
						$GENOME_ASSEMBLIES{$key} = $assembly;
					} elsif($path_type eq 'bowtie_path'){
						$BOWTIE_PATH_CONFIGS{$key} = $genome_file;
						$BOWTIE_PATHS{$key} = $path;
						$BOWTIE_SPECIES{$key} = $species;
						$BOWTIE_ASSEMBLIES{$key} = $assembly;
					} elsif($path_type eq 'gtf_path'){
						$GTF_PATH_CONFIGS{$key} = $genome_file;
						$GTF_PATHS{$key} = $path;
						$GTF_SPECIES{$key} = $species;
						$GTF_ASSEMBLIES{$key} = $assembly;
					}
				}
			}
			close(GCONFIG);
		}
	}
}


sub parse_updates_file {
	my $updates_file = $ENV{"HOME"}."/clusterflow/.cfupdates";
	if(-e $updates_file){
		open (UPDATES, $updates_file) or die "Can't read $updates_file: $!";
		$AVAILABLE_VERSION = <UPDATES>;
		$AVAILABLE_VERSION =~ s/[\n\r]//;
		$UPDATES_LAST_CHECKED = <UPDATES>;
		$UPDATES_LAST_CHECKED =~ s/[\n\r]//;
		close (UPDATES);
	}
}


# Lists available genomes
sub list_clusterflow_genomes {
	
	my $returnstring = "";
	
	my @config_files = ("$FindBin::Bin/genomes.config", "$homedir/clusterflow/genomes.config", './genomes.config');
	
	# See if we have any paths in each config file
	my %GENOME_PATH_CONFIGS_COUNTS;
	while ( my($key, $value) = each %GENOME_PATH_CONFIGS){
		$GENOME_PATH_CONFIGS_COUNTS{$value} = 1;
	}
	my %BOWTIE_PATH_CONFIGS_COUNTS;
	while ( my($key, $value) = each %BOWTIE_PATH_CONFIGS){
		$BOWTIE_PATH_CONFIGS_COUNTS{$value} = 1;
	}
	my %GTF_PATH_CONFIGS_COUNTS;
	while ( my($key, $value) = each %GTF_PATH_CONFIGS){
		$GTF_PATH_CONFIGS_COUNTS{$value} = 1;
	}
	
	foreach my $config_file (@config_files){
		
		next if(!defined ($GENOME_PATH_CONFIGS_COUNTS{$config_file}) && !defined ($BOWTIE_PATH_CONFIGS_COUNTS{$config_file}) && !defined ($GTF_PATH_CONFIGS_COUNTS{$config_file}));
		
		$returnstring .= "\n".('-' x 50)."\n $config_file\n".('-' x 50)."\n";
		if(defined ($GENOME_PATH_CONFIGS_COUNTS{$config_file})){
			$returnstring .= "\n== Genome Paths ==\n";
			$returnstring .= " Key                 Species             Assembly            Path\n".("-" x 100)."\n";
			foreach my $key (sort keys %GENOME_PATHS ) {
				my $key_spaces = " " x (20 - length($key));
				my $species_spaces = " " x (20 - length($GENOME_SPECIES{$key}));
				my $assembly_spaces = " " x (15 - length($GENOME_ASSEMBLIES{$key}));
				if($GENOME_PATH_CONFIGS{$key} eq $config_file){
					$returnstring .= " ".$key.$key_spaces.$GENOME_SPECIES{$key}.$species_spaces.$GENOME_ASSEMBLIES{$key}.$assembly_spaces.$GENOME_PATHS{$key}."\n";
				}
			}
		}
		if(defined ($BOWTIE_PATH_CONFIGS_COUNTS{$config_file})){
			$returnstring .= "\n== Bowtie Index Base Paths ==\n";
			$returnstring .= " Key                 Species             Assembly            Path\n".("-" x 100)."\n";
			foreach my $key (sort keys %BOWTIE_PATHS ) {
				my $key_spaces = " " x (20 - length($key));
				my $species_spaces = " " x (20 - length($BOWTIE_SPECIES{$key}));
				my $assembly_spaces = " " x (15 - length($BOWTIE_ASSEMBLIES{$key}));
				if($BOWTIE_PATH_CONFIGS{$key} eq $config_file){
					$returnstring .= " ".$key.$key_spaces.$BOWTIE_SPECIES{$key}.$species_spaces.$BOWTIE_ASSEMBLIES{$key}.$assembly_spaces.$BOWTIE_PATHS{$key}."\n";
				}
			}
		}
		if(defined ($GTF_PATH_CONFIGS_COUNTS{$config_file})){
			$returnstring .= "\n== GTF File Paths ==\n";
			$returnstring .= " Key                 Species             Assembly            Path\n".("-" x 100)."\n";
			foreach my $key (sort keys %GTF_PATHS ) {
				my $key_spaces = " " x (20 - length($key));
				my $species_spaces = " " x (20 - length($GTF_SPECIES{$key}));
				my $assembly_spaces = " " x (15 - length($GTF_ASSEMBLIES{$key}));
				if($GTF_PATH_CONFIGS{$key} eq $config_file){
					$returnstring .= " ".$key.$key_spaces.$GTF_SPECIES{$key}.$species_spaces.$GTF_ASSEMBLIES{$key}.$assembly_spaces.$GTF_PATHS{$key}."\n";
				}
			}
		}
	}
	$returnstring .= "\n";
	
	return $returnstring;
}


# Prints help for a specific module or pipeline
sub clusterflow_pipeline_help {
	
	my ($pipeline) = @_;
	
	my $help = "";
	
	my @pipelines = ("./$pipeline.config", "$homedir/clusterflow/pipelines/$pipeline.config", "$FindBin::Bin/pipelines/$pipeline.config");
	my @modules = ("./$pipeline.cfmod", "$homedir/clusterflow/modules/$pipeline.cfmod", "$FindBin::Bin/modules/$pipeline.cfmod");
	foreach my $pipeline (@pipelines){
		if(-e $pipeline){
			open (PIPELINE, $pipeline) or die "Can't read $pipeline: $!";
			my $comment_block = 0;
			while (<PIPELINE>) {
				chomp;
				s/\n//;
				s/\r//;
				next if($_ =~ /^\/\*/); # multiline comments start
				if($_ =~ /^\*\//){		# multiline comments end
					$help .= "\n".("-" x 20)."\n Pipeline:\n".("-" x 20)."\n";
					next;
				}
				$help .= $_."\n";
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
	
	Note that the name of a single module can be used instead of a
	pipeline name.

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
	
	--no_fn_check
		Disable input file type checking
	
	--file_list
		Text file containing input files or download URLs
		
	--mod_params
		Specify module parameters when running a single module.
		Ignored if a pipeline name is specified.
		
	--split_files <num>
		Create one run per <num> files
	
	--max_runs <num>
		Divide input files into <num> runs.
		Overrides --split_files
		
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
		
	--dry_run
		Prints jobs to terminal instead of submitting them to the cluster
	
	--qstat
		Parses the output from qstat in a visually attractive and intuitive manner
		
	--qstatall
		Same as --qstat, but for all jobs submitted by all users
	
	--qstatcols
		Colour the output from --qstat and --qstatall. Looks nice in terminal
		windows with light backgrounds. Horrible on dark backgrounds.
		See the manual for recommended .bashrc aliases.
	
	--qdel
		Delete all jobs running in a particular Cluster Flow pipeline. Follow
		with a pipeline ID, printed when running --qstat
		
	--make_config
		Interactive prompt to generate a personalised CF config file
	
	--version
		Print version of Cluster Flow installed
	
	--check_updates
		Look for available Cluster Flow updates
		
	--help
		Print this help message.
		If specified with a pipeline or module name afterwards, the help for that
		pipeline or module will be displayed. eg. cf --help sra_bismark


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
	
	my $use_defaults;
	my $use_defaults_stdin;
	print "Ok, the rest of this wizard is about which notification e-mails that
you'd like to receive. We can skip this and use default settings if you
prefer. Use defaults? (y/n)\n\n";
	while ($use_defaults_stdin = <STDIN>){
		chomp ($use_defaults_stdin);
		if($use_defaults_stdin =~ /^y(es)?/i){
			print "\nGood choice. You can always edit these later anyway, just see the manual..\n\n";
			$use_defaults = 1;
			sleep(2);
			last;
		} elsif ($use_defaults_stdin =~ /^n(o)?/i){
			print "\nOk, let's delve a little deeper..\n\n";
			last;
		} else {
			print "\nSorry, I didn't understand that.\nCould you try again please? (y/n)\n\n";
		}
	}
	if($use_defaults){
		$config .= "\@notification	complete\n";
		$config .= "\@notification	suspend\n";
		$config .= "\@notification	abort\n";
	} else {
		
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
	
	} # end of defaults check
	
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