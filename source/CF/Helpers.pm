#!/usr/bin/env perl
package CF::Helpers;

use warnings;
use strict;
use Exporter;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptionsFromArray);
use POSIX qw(ceil strftime);
use Time::Local;
use CF::Constants;

##########################################################################
# Copyright 2014, Philip Ewels (phil.ewels@scilifelab.se)                #
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


# Used by modules when called. Either prints the help text, the module
# requirements or returns a hash with information for the run.
# Takes three variables - the @ARGV reference, a hash with keys for
# cores, memory, modules and time requirements (can be arrays, strings or
# subroutines) and the help text string.
sub module_start {
    # Incoming..
    my ($clargs, $reqs, $helptext) = @_;

    # Get Command Line Options
    my $get_requirements;
    my @run_fns;
    my $job_id;
    my $prev_job_id;
    my $cores;
    my $mem;
    my %params;
    my $help;
    my $result = GetOptionsFromArray ($clargs,
        "requirements"      => \$get_requirements,
        "run_fn=s"          => \@run_fns,
        "job_id=s"          => \$job_id,
        "prev_job_id=s"     => \$prev_job_id,
        "cores=i"           => \$cores,
        "mem=s"             => \$mem,
        "param=s"           => \%params,
        "help"              => \$help
    );

    # For non-summary modules, there will only be one run file
    my $run_fn = $run_fns[0];

    # Get caller details
    my ($package, $mod_fn, $line) = caller;
    my $modname;
    ($modname = $mod_fn) =~ s/\.cfmod.pl$//i; # Always a perl file here
    $modname =~ s/^.*\///i;

    # Initialise the run file hash
    my %runfile = (
        'run_fn'        => $run_fn,
        'run_fns'       => \@run_fns,
        'modname'       => $modname,
        'mod_fn'        => $mod_fn,
        'job_id'        => $job_id,
        'prev_job_id'   => $prev_job_id,
        'cores'         => $cores,
        'memory'        => $mem,
        'params'        => \%params,
    );

    # Print help
    if($help){
        print $helptext;
        exit;
    }

    # Check that we have a run file, need it to go any further
    if(!$run_fn or length($run_fn) == 0){
        die ("Error: No run file filename supplied for module $modname\n");
    }

    ###
    # Module requirements
    ###
    if($get_requirements){

        # Parse the run file(s)
        parse_runfile(\%runfile);

        # Run through the supplied requirements and print
        foreach my $key (keys (%{$reqs})){
            # Been given a function
            if(ref($reqs->{$key}) eq 'CODE'){
                print "$key: ".$reqs->{$key}(\%runfile)."\n";
            }
            # Been given an array
            elsif(ref($reqs->{$key}) eq 'ARRAY'){
                if($key eq 'cores'){
                    print "cores: ".allocate_cores($cores, int($reqs->{$key}[0]), int($reqs->{$key}[1]))."\n";
                } elsif($key eq 'memory'){
                    print "memory: ".bytes_to_human_readable(allocate_memory($mem, $reqs->{$key}[0], $reqs->{$key}[1]))."\n";
                } else {
                    print "$key: ".join(', ', @{$reqs->{$key}})."\n";
                }
            }
            # Probably been given a string
            else{
                print "$key: ".$reqs->{$key}."\n";
            }
        }
        exit;
    }

    ###
    # Running for real
    ###
    # check that we have extra parameters that we need
    if(!$job_id or length($job_id) == 0 or !$prev_job_id or length($prev_job_id) == 0){
        die ("###CF Error: Need job ID and previous job ID! Missing for module $modname\n");
    }
    if(!$cores or length($cores) == 0 or !$mem or length($mem) == 0){
        die ("###CF Error: Need allocated cores and memory! Missing for module $modname\n");
    }

    # Parse everything we can from the run file
    parse_runfile(\%runfile);

    # Print the module header if needed
    if(!defined($params{'hide_log_header'})){
        my $summ = '';
        $summ = "Summary Module:\t\tYes\n" if(defined($params{'summary_module'}));

        my $params = '';
        while( my( $key, $value ) = each %params ){
            $params .= "$key = $value, ";
        }
        $params = "Parameters:\t\t".substr($params, 0, -2)."\n" if(length($params) > 0);

    	my $date = strftime "%H:%M, %d-%m-%Y", localtime;
    	my $dashes = "-" x 80;

    	warn "\n$dashes\nModule:\t\t\t$modname\n".$summ."Run File:\t\t$run_fn\n";
        warn "Job ID:\t\t\t$job_id\nPrevious Job ID:\t$prev_job_id\n";
        warn $params."Date & Time:\t\t$date\n$dashes\n\n";
    }

	# If we don't have any input files, bail now
	if(!defined($runfile{'prev_job_files'}) || scalar(@{$runfile{'prev_job_files'}}) == 0){
        if($prev_job_id && $prev_job_id ne 'null' && !defined($params{'summary_module'})){
    	    print "\n###CF Error! No file names found from job $prev_job_id. Exiting...\n\n";
    	    exit;
        }
	}

    # Return the hash of dreams
    return %runfile;

}


# Open the run file and parse it's contents.
# Input arguments: Reference to %runfile
# Returns: Nothing (updates hash reference in place)
sub parse_runfile {

    # Get incoming hash reference
    my $runfile = $_[0];

    # Set up new hash variables if we need them
    $runfile->{'refs'} = {}                             if(!defined($runfile->{'refs'}));
	$runfile->{'config'} = {}                           if(!defined($runfile->{'config'}));
	$runfile->{'config'}{'notifications'} = {}          if(!defined($runfile->{'config'}{'notifications'}));
	$runfile->{'prev_job_files'} = ()                   if(!defined($runfile->{'prev_job_files'}));
	$runfile->{'starting_files'} = ()                   if(!defined($runfile->{'starting_files'}));
    $runfile->{'files'} = {}                            if(!defined($runfile->{'files'}));
    $runfile->{'num_starting_files'} = 0                if(!defined($runfile->{'num_starting_files'}));
    $runfile->{'num_starting_merged_files'} = 0         if(!defined($runfile->{'num_starting_merged_files'}));
    $runfile->{'num_starting_merged_aligned_files'} = 0 if(!defined($runfile->{'num_starting_merged_aligned_files'}));

    # Go through each run file
    foreach my $runfn (@{$runfile->{'run_fns'}}){
        open (RUN, $runfn) or die "Can't read $runfn: $!";
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

            # Helper stuff
            if($_ =~ /^Pipeline: (.+)$/){
                $runfile->{'pipeline_name'} = $1;
            }

    		# Get config variables
    		if($_ =~ /^\@/ && !$comment_block){
    			my @sections = split(/\t/, $_, 2);
    			my $cname = substr($sections[0], 1);
                $sections[1] = 1 if(!defined($sections[1]));
    			if($cname eq 'notification'){
    				$runfile->{'config'}{'notifications'}{$sections[1]} = 1;
    			} elsif($cname eq 'reference'){
                    my @ref_sections = split(/\t/, $sections[1], 2);
    				$runfile->{'refs'}{$ref_sections[0]} = $ref_sections[1];
    			} else {
    				$runfile->{'config'}{$cname} = $sections[1];
    			}
    		}

            # Note - could get the pipeline tree here. Currently no need.

    		# Get files
    		if($_ =~ /^[^@#]/ && !$comment_block){
    			my @sections = split(/\t/, $_, 2);

                # Clear out excess whitespace
                $sections[1] =~ s/^\s+//;
                $sections[1] =~ s/\s+$//;

                # Files from previous job
    			if(defined($runfile->{'prev_job_id'}) && $sections[0] eq $runfile->{'prev_job_id'}){
    				push(@{$runfile->{'prev_job_files'}}, $sections[1]);
    			}

                # Starting files
                if($sections[0] eq 'start_000'){
                    push(@{$runfile->{'starting_files'}}, $sections[1]);
                    $runfile->{'num_starting_files'}++;
                }

                # All files, by module
                if(!defined($runfile->{'files'}{$sections[0]})){
                    $runfile->{'files'}{$sections[0]} = ();
                }
                push(@{$runfile->{'files'}{$sections[0]}}, $sections[1]);

    		}
    	}

    	close(RUN);
    }

    # Figure out how many merged files we're likely to have
    my $regex;
    if(defined($runfile->{'config'}{'merge_regex'}) && length($runfile->{'config'}{'merge_regex'}) > 0){
    	$regex = $runfile->{'config'}{'merge_regex'};
    }
    if(defined($runfile->{'params'}{'regex'})){
        $regex = $runfile->{'params'}{'regex'};
    }
    if($regex){
        my %file_sets;
        for my $file (@{$runfile->{'starting_files'}}){
        	my $group = ($file =~ m/$regex/) ? $1 : 'file';
        	if(!defined($file_sets{$group})){
        		$file_sets{$group} = 1;
        	}
        }
        $runfile->{'num_starting_merged_files'} = scalar(keys(%file_sets));
    } else {
        $runfile->{'num_starting_merged_files'} = $runfile->{'num_starting_files'};
    }

    # How many aligned files? Just check for paired end or single end config
    if(defined($runfile->{'config'}{'force_paired_end'})){
        $runfile->{'num_starting_merged_aligned_files'} = $runfile->{'num_starting_merged_files'} / 2;
    } else {
        $runfile->{'num_starting_merged_aligned_files'} = $runfile->{'num_starting_merged_files'};
    }

    # No need to return anything, as we've been using a hash reference
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
			my $mod_cmd = `modulecmd perl load $mod 2> /dev/null`;
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

	my $runfile = shift;

	my @files = sort(@_);
	my @se_files;
	my @pe_files;

	# Force Paired End or Single End if specified in the config
	if(exists($runfile->{'config'}{'force_paired_end'})){
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
	} elsif(exists($runfile->{'config'}{'force_paired_end'})){
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
# Reads through the first 1000 reads and decides on how many 0x1 flags it finds
sub is_bam_paired_end {

	# Load samtools
	my @modules = ('samtools');
	my %loaded_mods = ();
	&CF::Helpers::load_environment_modules(\@modules,\%loaded_mods);

	my ($file) = @_;

	unless($file =~ /.bam$/ || $file =~ /.sam$/){
		warn "\n$file is not a .bam or .sam file - can't figure out PE / SE mode..\nExiting..\n\n";
		die;
	}

    # Read the first 1000 lines withx samtools
    my $se_reads = 0;
    my $pe_reads = 0;
    my $readcount = 0;
    open(my $fh, '-|', "samtools view $file") or die "Could not run samtools to check BAM pe/se: $!";
    while (<$fh>) {
        last if ($readcount >= 1000);
        my ($flag) = (split (/\t/))[1];
        if ($flag & 0x1){
            $pe_reads++;
        } else {
            $se_reads++;
        }
        $readcount++;
    }
    close($fh);

    # Look at our counts
    if($pe_reads >= 800){
        return 1;
    } else {
        return 0;
    }

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


# Function to convert a SLURM style timestamp to minutes
# minutes
# minutes:seconds
# hours:minutes:seconds
# days-hours
# days-hours:minutes
# days-hours:minutes:seconds
sub timestamp_to_minutes {
    my $timestamp = $_[0];
    my $days = 0;
    my $hours = 0;
    my $minutes = 0;
    my $seconds = 0;
    my $re_hh = "([0-2]?[0-9])";
    my $re_mm = "([0-5]?[0-9])";

    # Make sure that we don't have any illegal characters
    if($timestamp =~ /[^\d\-\:]/){
        return 0;
    }

    # First - days
    if($timestamp =~ /^(\d+)\-/){
        $days = $1;
        $timestamp =~ s/^(\d+)\-//;
    }

    # Progressively match colons to figure out what the rest is
    if($timestamp =~ /^$re_hh\:$re_mm\:$re_mm$/){
        $hours = $1;
        $minutes = $2;
        $seconds = $3;
    } else {
        if($days > 0){
            if($timestamp =~ /^$re_hh\:$re_mm$/){
                $hours = $1;
                $minutes = $2;
            } elsif(/^$re_hh$/){
                $hours = $1;
            }
        } else {
            if($timestamp =~ /^$re_mm\:$re_mm$/){
                $minutes = $1;
                $seconds = $2;
            } elsif($timestamp =~ /^$re_mm$/){
                $minutes = $1;
            }
        }
    }

    my $total_minutes = 0;
    $total_minutes += $days * 24 * 60;
    $total_minutes += $hours * 60;
    $total_minutes += $minutes;
    # We ignore seconds. Seriously, who does that?

    return $total_minutes;

}

# Function to convert minutes to a SLURM time stamp. Reverse of above.
sub minutes_to_timestamp {
    my $minutes = $_[0];

    # Check for illegal characters
    if($minutes =~ /[^\d\.]/){
        return 0;
    }
    $minutes = int($minutes);

    my $days = int($minutes / (24 * 60));
    $minutes -= $days * 24 * 60;
    $minutes = ($minutes < 0) ? 0 : $minutes;

    my $hours = int($minutes / 60);
    $minutes -= $hours * 60;
    $minutes = ($minutes < 0) ? 0 : $minutes;

    $minutes = sprintf("%02d", $minutes);

    my $timestamp = '';
    if($days > 0){
        return "$days-$hours:$minutes";
    } elsif($hours > 0){
        return "$hours:$minutes:00";
    } else {
        return "$minutes";
    }
}





# Function to take human readable memory string and return bytes
sub human_readable_to_bytes {

	my $memory = $_[0];
    my $suffix = '';
    if($memory =~ /([tgmk])b?$/i){
        $suffix = $1;
    }
	$memory =~ s/[^\d\.]//g;

	if(lc($suffix) eq 't'){
		$memory = $memory * 1000000000000;
	} elsif(lc($suffix) eq 'g'){
		$memory = $memory * 1000000000;
	} elsif(lc($suffix) eq 'm'){
		$memory = $memory * 1000000;
	} elsif(lc($suffix) eq 'k'){
		$memory = $memory * 1000;
	}

	return $memory;
}

# Function to take bytes and return a human readable memory string
sub bytes_to_human_readable {

	my ($bytes) = @_;
	$bytes =~ s/\D//g;

	if(int($bytes/1000000000) > 0){
		return ceil($bytes/1000000000)."G";
	} elsif(int($bytes/1000000) > 0){
		return ceil($bytes/1000000)."M";
	} elsif(int($bytes/1000) > 0){
		return ceil($bytes/1000)."K";
	} else {
		return $bytes."B";
	}

}

# Simple function to take bytes and return a human readable memory string
# NB: Rounds up with ceil()
sub mem_return_mbs {

	my ($mem) = @_;
	$mem = human_readable_to_bytes($mem);
	return ceil($mem/1000000);

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

    $allocated = human_readable_to_bytes($allocated);
	$max = human_readable_to_bytes($max);
	$min = human_readable_to_bytes($min);

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
sub build_emails {

  my ($title, $html_content, $plain_content) = @_;

	# Get the e-mail templates
	# Assume that we're running from the installation directory/modules
	my $html_email;
	{ local $/ = undef; local *FILE; open FILE, "<$Bin/../source/CF/html_email_template.html"; $html_email = <FILE>; close FILE }
	my $text_email;
	{ local $/ = undef; local *FILE; open FILE, "<$Bin/../source/CF/plaintext_email_template.txt"; $text_email = <FILE>; close FILE }

	my $cf_version = $CF::Constants::CF_VERSION;

	# Put in our content
	$html_email =~ s/{{ PAGE_TITLE }}/$title/g;
	$html_email =~ s/{{ CONTENT }}/$html_content/g;
	$html_email =~ s/{{ CF_VERSION }}/$cf_version/g;

	$text_email =~ s/{{ PAGE_TITLE }}/$title/g;
	$text_email =~ s/{{ CONTENT }}/$plain_content/g;
	$text_email =~ s/{{ CF_VERSION }}/$cf_version/g;

  return($html_email, $text_email);

}

sub send_email {

	my ($subject, $to, $title, $html_content, $plain_content) = @_;

  my ($html_email, $text_email) = build_emails($title, $html_content, $plain_content);

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
