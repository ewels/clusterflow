#!/usr/bin/env perl
package CF::Helpers; 

use warnings;
use strict;
use FindBin qw($Bin);
use Exporter;
use POSIX qw(strftime);
use Time::Local;
use CF::Constants; 

##########################################################################
# Copyright 2014, Philip Ewels (phil.ewels@babraham.ac.uk)               #
#                                                                        #
# This file is part of Cluster Flow.                                     #
#                                                                        #
# Cluster Flow is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# Cluster Flow is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with Cluster Flow.  If not, see <http://www.gnu.org/licenses/>.  #
##########################################################################


sub load_runfile_params {
	my ($runfile, $job_id, $prev_job_id, $cores, $mem, @parameters) = @_;
	unless (defined $prev_job_id && length($prev_job_id) > 0) {
		die "Previous job ID not specified\n";
	}
	unless (@parameters) {
		@parameters = ();
	}
	
	my $num_input_files = 0;
	
	my $date = strftime "%H:%M, %d-%m-%Y", localtime;
	my $dashes = "-" x 80;
	warn "\n$dashes\nRun File:\t\t$runfile\nJob ID:\t\t\t$job_id\nPrevious Job ID:\t$prev_job_id\nParameters:\t\t".join(", ", @parameters)."\nDate & Time:\t\t$date\n$dashes\n\n";

	open (RUN,$runfile) or die "Can't read $runfile: $!";

	my @files;
	my %config;
	$config{notifications} = {};
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
			my @sections = split(/\t/, $_, 2);
			my $cname = substr($sections[0], 1);
			if($cname eq 'notification'){
				$config{notifications}{$sections[1]} = 1;
			} else {
				$config{$cname} = $sections[1];
			}
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
	if($num_input_files == 0 && $prev_job_id ne 'null'){
		print "\n###CF Error! No file names found from job $prev_job_id. Exiting...\n\n";
		exit;
	}
	
	return (\@files, $runfile, $job_id, $prev_job_id, $cores, $mem, \@parameters, \%config);
}


sub parse_runfile_prerun {
	my $runfile = @_[0];
	my $num_input_files = 0;
	my @starting_files;
	my %config;
	$config{notifications} = {};
	my $comment_block = 0;
	
	open (RUN,$runfile) or die "Can't read $runfile: $!";
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
			my @sections = split(/\t/, $_, 2);
			my $cname = substr($sections[0], 1);
			if($cname eq 'notification'){
				$config{notifications}{$sections[1]} = 1;
			} else {
				$config{$cname} = $sections[1];
			}
		}
		
		# Get files
		if($_ =~ /^[^@#]/ && !$comment_block){
			my @sections = split(/\t/, $_, 2);
			if($sections[0] eq "start_000"){
				# Clear out excess whitespace
				$sections[1] =~ s/^\s+//;
				$sections[1] =~ s/\s+$//;
				# Push to array
				push(@starting_files, $sections[1]);
				$num_input_files++;
			}
		}
	}
	
	close(RUN);
	
	# If we don't have any input files, bail now
	if($num_input_files == 0){
		warn "\n###CF Error! No starting file names found (runfile should contain start_000  filename). Exiting...\n\n";
		exit;
	}
	return(\@starting_files, \%config)
}




# Function to load environment modules into the perl environment
sub load_environment_modules {
	my ($modules, $loaded_modules) = @_;
	my $use_modules = $CF::Constants::CF_MODULES;
	my %mod_aliases = %CF::Constants::ENV_MODULE_ALIASES;
	if($CF::Constants::CF_MODULES){
		foreach my $mod (@{$modules}) {
			# Check to see if we have an alias for this module
			if(defined($mod_aliases{$mod})){
				$mod = $mod_aliases{$mod};
			}
			# Skip modules that have already been loaded
			next if defined($loaded_modules->{$mod});
			# Get perl code needed to load module
			my $mod_cmd = `modulecmd perl load $mod`;
			if($mod_cmd && length($mod_cmd) > 0){
				eval($mod_cmd);
				if ($@){
					warn "WARNING - Got error whilst trying to parse the module load code for module $mod:" .
					"\t$@\n\nTHIS MODULE HAS NOT BEEN LOADED. Skipping..\n\n";
					sleep(2);
				} else {
					# Everything worked. Remember so we don't try again.
					$loaded_modules->{$mod} = 1;
				}
			}
		}
	}
	return (\%{$loaded_modules});
}




# Function to look at supplied file names and work out whether they're paired end or not
sub is_paired_end {
	
	my $config = shift;
	
	my @files = sort(@_);
	my @se_files;
	my @pe_files;
	
	# Force Paired End or Single End if specified in the config
	if(exists($config->{force_paired_end})){
		for (my $i = 0; $i <= $#files; $i++){
			if($i < $#files){
				my @pe = ($files[$i], $files[$i+1]);
				push (@pe_files, \@pe);
				$i++;
			} else {
				# If we have an odd number of file names and have got to the end, must be SE
				push (@se_files, $files[$i]);
			}
		}
		return (\@se_files, \@pe_files);
	} elsif(exists($config->{force_single_end})){
		for (my $i = 0; $i <= $#files; $i++){
			push (@se_files, $files[$i]);
		}
		return (\@se_files, \@pe_files);
	}
	
	# Haven't returned yet, so let's figure it out for ourselves
	for (my $i = 0; $i <= $#files; $i++){
		if($i < $#files){
			# Make stripped copies of the fns for comparison
			(my $fn1 = $files[$i]) =~ s/_R?[1-4]//g;
			(my $fn2 = $files[$i+1]) =~ s/_R?[1-4]//g;
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


# Function to look into BAM/SAM header to see whether it's paired end or not
sub is_bam_paired_end {

	# Load samtools
	my @modules = ('samtools');
	my %loaded_mods = {};
	&CF::Helpers::load_environment_modules(\@modules,\%loaded_mods);

	my ($file) = @_;
	
	unless($file =~ /.bam$/ || $file =~ /.sam$/){
		warn "\n$file is not a .bam or .sam file - can't figure out PE / SE mode..\nExiting..\n\n";
		die;
	}
	
	my $headers = `samtools view -H $file`;
	
	my $paired_end;
	my $header_found;
	foreach (split("\n", $headers)){
		if(/^\@PG/){
			if(/-1/ && /-2/){
				$paired_end = 1;
				$header_found = 1;
			} else {
				$paired_end = 0;
				$header_found = 1;
			}
		}
	}
	
	if(!$header_found){
		warn "\nCould not find BAM \@PG header, so could not determine read type.\nExiting...\n\n";
		die;
	}
	
	return $paired_end;
	
}

# Function to determine the encoding of a FastQ file (nicked from HiCUP)
# See http://en.wikipedia.org/wiki/FASTQ_format#Encoding
# phred33: ASCII chars begin at 33
# phred64: ASCII chars begin at 64
# solexa: ASCII chars begin at 59
# integer: quality values integers separated by spaces
sub fastq_encoding_type {

	my $file = $_[0];
	my $score_min = 999;    #Initialise at off-the-scale values
	my $score_max = -999;
	my $read_count = 0;
	
	if($file =~ /\.gz$/){
		open (IN, "zcat $file |") or die "Could not read file '$file' : $!";
	} else {
		open (IN, $file) or die "Could not read file '$file' : $!";
	}
	
	while(<IN>){
	
		unless(/^@/){
			# Line must start with an @ symbol - read identifiers
			die "Error trying to work out the FastQ quality scores!\nRead doesn't start with an \@ symbol..\n\n\n";
		}
			
		# push file counter on two lines to the quality score
		scalar <IN>;
		scalar <IN>; 
		
		my $quality_line = scalar <IN>;
		chomp $quality_line;    
		my @scores = split(//, $quality_line);
		
		foreach(@scores){
			my $score = ord $_;    #Determine the value of the ASCII character
			
			if($score < $score_min){
				$score_min = $score;
			}
			if($score > $score_max){
				$score_max = $score;
			}	
		}
		
		# Do not need to process 100,000 lines if these parameters are met
		if($score_min == 32){    # Contains the space charcter
			close IN;
			return 'integer';
		} elsif ($score_min < 59){   # Contains low range character
			close IN;
			return 'phred33';
		} elsif ( ($score_min < 64) and ($score_max > 90) ){	# Contains character below phred64 and far above phred33
			close IN;
			return 'solexa'
		}
		
		$read_count++;
		
	}
	close IN;
	
	if($read_count < 100000){
		return 0;    #File did not contain enough lines to make a decision on quality
	} else {
		return 'phred64';
	}
}


# Function to determine the maximum read length
# found in the first 100000 lines of a FastQ file
sub fastq_min_length {

	my ($file, $minlength) = @_;
	my $read_count = 0;
	
	if($file =~ /\.gz$/){
		open (IN, "zcat $file |") or die "Could not read file '$file' : $!";
	} else {
		open (IN, $file) or die "Could not read file '$file' : $!";
	}
	
	while(<IN>){
		
		unless(/^@/){
			# Line must start with an @ symbol - read identifiers
			die "Error trying to work out the FastQ read lengths!\nRead doesn't start with an \@ symbol..\n\n\n";
		}
		
		# push file counter on two lines to the quality score
		scalar <IN>;
		scalar <IN>;
		
		my $quality_line = scalar <IN>;
		chomp $quality_line;
		
		if(length($quality_line) >= $minlength){
			close IN;
			return 1;
		}
		
		$read_count++;
		
		if($read_count > 100000){
			last;
		}
		
	}
	close IN;
	
	return 0;
	
}


# Simple function to take time in seconds and convert to human readable string
sub parse_seconds {

	my ($raw, $long) = @_;
	unless(defined($long)){
		$long = 1;
	}
	
	
	my @chunks;
	
	my $w_secs = 's';
	my $w_mins = 'm';
	my $w_hours = 'h';
	my $w_days = 'd';
	if($long){
		$w_secs = ' seconds';
		$w_mins = ' minutes';
		$w_hours = ' hours';
		$w_days = ' days';
	}
	
	my $days = int($raw/(24*60*60));
	if($days > 0){
		push (@chunks, "$days$w_days");
		$raw -= $days * (24*60*60);
	}
	
	my $hours = ($raw/(60*60))%24;
	if($hours > 0){
		push (@chunks, "$hours$w_hours");
		$raw -= $hours * (60*60);
	}
	
	my $mins = ($raw/60)%60;
	if($mins > 0){
		push (@chunks, "$mins$w_mins");
		$raw -= $mins * 60;
	}
	
	my $secs = $raw%60;
	if($secs > 0){
		push (@chunks, "$secs$w_secs");
	}
	
	my $output = join(" ", @chunks);
	if($long){
		$output = join(", ", @chunks);
	}
	
	
	return ($output);
}



# Simple function to take human readable memory string and return bytes
sub human_readable_to_bytes {

	my ($memory_string) = @_;
	(my $memory = $memory_string) =~ s/\D//g;
	
	# Cut off 'b' from Gb etc if accidentally present
	$memory_string =~ s/b$//i;
	
	if(substr(lc($memory_string), -1) eq 'g'){
		$memory = $memory * 1073741824;
	} elsif(substr(lc($memory_string), -1) eq 'm'){
		$memory = $memory * 1048576;
	} elsif(substr(lc($memory_string), -1) eq 'k'){
		$memory = $memory * 1024;
	}
	
	return $memory;
}

# Simple function to take bytes and return a human readable memory string
sub bytes_to_human_readable {

	my ($bytes) = @_;
	$bytes =~ s/\D//g;
	
	if(int($bytes/1073741824) > 0){
		return sprintf("%.1f", $bytes/1073741824)."G";
	} elsif(int($bytes/1048576) > 0){
		return sprintf("%.1f", $bytes/1048576)."M";
	} elsif(int($bytes/1024) > 0){
		return sprintf("%.1f", $bytes/1048576)."K";
	} else {
		return $bytes."B";
	}
	
}

# Simple function to take bytes and return a human readable memory string
sub mem_return_mbs {

	my ($mem) = @_;
	$mem = human_readable_to_bytes($mem);
	return sprintf("%.0f", $mem/1048576);
	
}

# Take allocated cores, minimum, maximum and return best value
sub allocate_cores {
	
	my ($allocated, $min, $max) = @_;
	
	if($allocated > $max){
		return $max;
	} elsif($allocated < $min){
		return $min;
	} else {
		return $allocated;
	}
}

# Take allocated memory, minimum, maximum and return best value
sub allocate_memory {
	
	my ($allocated, $min, $max) = @_;
	
	$max = CF::Helpers::human_readable_to_bytes($max);
	$min = CF::Helpers::human_readable_to_bytes($min);
	
	if($allocated > $max){
		return $max;
	} elsif($allocated < $min){
		return $min;
	} else {
		return $allocated;
	}
}



# Function to properly split up version numbers
# Returns true if second supplied vn is newer than first
sub cf_compare_version_numbers {
	
	my ($vn1_string, $vn2_string) = @_;
	
	my @vn1_parts = split(/\.|\s+/, $vn1_string);
	my @vn2_parts = split(/\.|\s+/, $vn2_string);
	
	for my $i (0 .. $#vn2_parts){
		if(defined($vn1_parts[$i])){
			
			# Numeric checks
			if($vn1_parts[$i] =~ /^\d+$/ && $vn2_parts[$i] =~ /^\d+$/){
				if($vn2_parts[$i] > $vn1_parts[$i]){
					return 1;
				}
			} elsif($vn1_parts[$i] !~ /^\d+$/ && $vn2_parts[$i] =~ /^\d+$/){
				# 0.1.1 beats 0.1 devel
				return 1;
			}
		} else {
			return 1;
		}
	}
	
	return 0;
	
}


### E-MAIL FUNCTIONS
sub send_email {
	
	my ($subject, $to, $title, $html_content, $plain_content) = @_;
	
	my $cf_version = $CF::Constants::CF_VERSION;
	
	# Get the e-mail templates
	# Assume that we're running from the installation directory/modules
	my $html_email;
	{ local $/ = undef; local *FILE; open FILE, "<$Bin/../source/CF/html_email_template.html"; $html_email = <FILE>; close FILE }
	my $text_email;
	{ local $/ = undef; local *FILE; open FILE, "<$Bin/../source/CF/plaintext_email_template.txt"; $text_email = <FILE>; close FILE }
	
	# Put in our content
	$html_email =~ s/{{ PAGE_TITLE }}/$title/g;
	$html_email =~ s/{{ CONTENT }}/$html_content/g;
	$html_email =~ s/{{ CF_VERSION }}/$cf_version/g;
	
	$text_email =~ s/{{ PAGE_TITLE }}/$title/g;
	$text_email =~ s/{{ CONTENT }}/$plain_content/g;
	$text_email =~ s/{{ CF_VERSION }}/$cf_version/g;
	
	# Do we have the Perl modules that we need?
	my $mail;
	my $mail_packages = eval "use Email::MIME::CreateHTML; use Email::Sender::Simple qw(sendmail); 1;";
	
	# Send a fancy HTML e-mail using perl packages
	if($mail_packages){
		warn "Sending e-mail using Perl packages..\n";
		my $email = Email::MIME->create_html(
			header => [
				# From => 'my@address',
				To => $to,
				Subject => '[CF] '.$subject
			],
			body => $html_email,
			text_body => $text_email
		);
		Email::Sender::Simple->send($email);
	
	# We don't have them, try with sendmail
	} elsif (!system('which sendmail > /dev/null 2>&1')) {
		warn "Sending HTML e-mail with sendmail..\n";
		$html_email = "To: $to\nSubject: [CF] $subject\nMime-Version: 1.0\nContent-Type: text/html\n\n".$html_email;
		open (PIPE , "| sendmail -t") or die "can't open pipe to sendmail: $!\n";
		print PIPE $html_email;
		close PIPE;
	
	# We don't have them, try HTML with mailx
	} elsif (!system('which mailx > /dev/null 2>&1')) {
		warn "Sending HTML e-mail with mailx..\n";
		open (PIPE , "| mail  -s '$(echo -e \"[CF] $subject\nContent-type: text/html;\")' $to") or die "can't open pipe to mailx: $!\n";
		print PIPE $html_email;
		close PIPE;
	
	# Fallback - use the basic mail with plaintext
	} else {
		warn "Sending e-mail using basic plain text mail..\n";
		open (PIPE , "| mail -s '[CF] $subject' $to") or die "can't open pipe to mail: $!\n";
		print PIPE $text_email;
		close PIPE;
	}
	
	# Give the program time to send the e-mail before qsub shuts us down
	sleep(5); 
	
	return 1;
	
}



1; # Must return a true value