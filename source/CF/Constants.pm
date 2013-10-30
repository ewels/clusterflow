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
	
	foreach (@NOTIFICATIONS){
		$output .= "\@notification\t$_\n";
	}
		
	return ($output);
	
}

# Prints help for a specific module or pipeline
sub clusterflow_pipeline_help {
	
	my ($pipeline) = @_;
	
	my $help;
	
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
		
	--p
		Force paired-end mode
		
	--s
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
		
	--notifications <eas>
		Specify desired notifications
		
	--list_pipelines
		Print available pipelines
		
	--list_modules
		Print available modules
		
	--list_genomes
		Print available genomes
		
	--dryrun
		Prints jobs to terminal instead of submitting them to the cluster
		
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


1;