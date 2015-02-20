#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;
use POSIX qw(strftime);
use Cwd;

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

# Get Options
my $required_cores;
my $required_mem;
my $required_modules;
my $run_fn;
my $help;
my $result = GetOptions ("cores=i" => \$required_cores, "mem=s" => \$required_mem, "modules" => \$required_modules, "runfn" => \$run_fn, "help" => \$help);
if($required_cores || $required_mem || $required_modules){
	exit;
}
if($help){
	die "\nThis is a core module which is executed when all runs have finished.\n";
}


# MODULE
# Read in the input files from the run file
my ($files, $runfile, $job_id, $prev_job_id, $cores, $mem, $parameters, $config_ref) = CF::Helpers::load_runfile_params(@ARGV);
my %config = %$config_ref;
@$parameters = grep { $_ ne 'summary_module' } @$parameters;
my $pipeline = shift(@$parameters);
my @outfns = @$parameters;

# Print run finish status to outfile
my $date = strftime ("%H:%M %d-%m-%Y", localtime);
warn "\n###CF Pipeline $pipeline finished at $date\n\n";

# Dig out pipeline ID and start from job ID
my $startdate = "?";
my $duration = "?";
my $pipeline_id = 'unknown_pipeline';
if($job_id =~ /^cf_(.+)_(\d{10})_(.+)_\d{1,3}$/){
	my $startdate_epoch = $2;
  $pipeline_id = "cf_$1_$2";
	$startdate = strftime ("%H:%M %d-%m-%Y", localtime($startdate_epoch));
	# Calculate duration
	my $duration_secs = time() - $startdate_epoch;
	$duration = CF::Helpers::parse_seconds($duration_secs);
}

# Find current directory
my $cwd = getcwd()."/";



# Read in the log files
my @LOG_HIGHLIGHT_STRINGS = @CF::Constants::LOG_HIGHLIGHT_STRINGS;
my @LOG_WARNING_STRINGS = @CF::Constants::LOG_WARNING_STRINGS;
my %cf_highlights;
my %commands;
my %warninglines;
my %highlightlines;
my @summarylines;
my $errors = 0;
my $warnings = 0;
my $highlights = 0;
foreach my $outfile (@outfns){
	
	my @these_cf_highlights;
	my @these_commands;
	my @these_highlightlines;
	my @these_warninglines;
    my @these_summarylines;
	
	open (IN,'<',$outfile);
	while(<IN>){
		
		chomp;
		
		# Ignore crap
		if(/^Warning: no access to tty/ || /^Thus no job control in this shell/){
			next;
		}
		
		# Commands run
		if(/^###CFCMD/){
			push (@these_commands, substr($_, 9));
		}
        
        # Summary statuses
        elsif(/^###CFSUMMARY/){
    	    push (@summarylines, substr($_, 13));
        }
        
		# Highlight statuses
		elsif(/^###CF/){
			push (@these_cf_highlights, substr($_, 6));
		} else {
			# Count any custom string highlights
			foreach my $highlight_string (@LOG_HIGHLIGHT_STRINGS){
				if(/$highlight_string/i){
					$highlights++;
					push (@these_highlightlines, $_);
				}
			}
		
			# Count any custom string errors
			foreach my $warning_string (@LOG_WARNING_STRINGS){
				if(/$warning_string/i){
					$warnings++;
					push (@these_warninglines, $_);
				}
			}
		}
        
		
		# Count out any CF errors
		if(/error/i){
			if (/^###CF/){
				$errors++;
			}
		}
		
	}
	close (IN);
	$cf_highlights{$outfile} = \@these_cf_highlights;
	$commands{$outfile} = \@these_commands;
	$highlightlines{$outfile} = \@these_highlightlines;
	$warninglines{$outfile} = \@these_warninglines;
}


# Send e-mail to submitter, if the config demands it
if($config{notifications}{complete} && defined($config{email})){
	
	# Write the plain-text e-mail body
my $plain_content = "The pipeline $pipeline has completed";
if($errors > 0){
	$plain_content .= " **with errors**";
} elsif($warnings > 0){
	$plain_content .= " **with warnings**";
}
$plain_content .= ".

Started: $startdate
Finished: $date
Duration: $duration
Working Directory: $cwd

";

if($warnings > 0){
	$plain_content .= "\n\n
===========================
== Custom Warnings Found ==
===========================
";
	foreach my $file (sort keys %warninglines) {
		if(scalar @{$warninglines{$file}} > 0){
			$plain_content .= "\n".("=" x 30)."\n- Output file $file\n".("=" x 30)."\n";
			$plain_content .= join("\n - ", @{$warninglines{$file}});
		}
	}	
}

if($highlights > 0){
	$plain_content .= "\n\n\n\n\n
=============================
== Custom Highlights Found ==
=============================
";
	foreach my $file (sort keys %highlightlines) {
		if(scalar @{$highlightlines{$file}} > 0){
			$plain_content .= "\n".("=" x 30)."\n- Output file $file\n".("=" x 30)."\n";
			$plain_content .= join("\n - ", @{$highlightlines{$file}});
		}
	}
	$plain_content .= "\n\n\n\n";
}

$plain_content .= "\n\n
=================================
== CF Status Messages ==
=================================
";
foreach my $file (sort keys %cf_highlights) {
	$plain_content .= "\n".("=" x 30)."\n- Output file $file\n".("=" x 30)."\n";
	$plain_content .= join("\n - ", @{$cf_highlights{$file}});
}
if(scalar @summarylines > 0){
	$plain_content .= "\n".("=" x 30)."\n- Summary Modules\n".("=" x 30)."\n";
	$plain_content .= join("\n - ", @summarylines);
}
$plain_content .= "\n\n\n\n\n
==================
== Commands Run ==
==================
";
foreach my $file (sort keys %commands) {
	$plain_content .= "\n".("=" x 30)."\n- Output file $file\n".("=" x 30)."\n";
	$plain_content .= join("\n\n", @{$commands{$file}});
}






	# Write the html e-mail body
	# Inline styles make me want to stab my eyes out, but we have to do this for
	# e-mail readers such as Gmail which strip any header CSS.
	my $html_content = '';
	if($errors > 0){
		$html_content .= '<p class="completion-leader" style="color: #723736; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: center; line-height: 19px; font-size: 18px; border-radius: 5px; background: #f2dede; margin: 0 0 15px; padding: 10px; border: 1px solid #ebccd1;" align="center">
The pipeline <strong>'.$pipeline.'</strong> has completed <strong>with errors</strong></p>';
	} elsif($warnings > 0){
		$html_content .= '<p class="completion-leader" style="color: #4D3E25; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: center; line-height: 19px; font-size: 18px; border-radius: 5px; background: #fcf8e3; margin: 0 0 15px; padding: 10px; border: 1px solid #faebcc;" align="center">
The pipeline <strong>'.$pipeline.'</strong> has completed <strong>with warnings</strong></p>';
	} else {
		$html_content .= '<p class="completion-leader" style="color: #0A440B; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: center; line-height: 19px; font-size: 18px; border-radius: 5px; background: #dff0d8; margin: 0 0 15px; padding: 10px; border: 1px solid #d6e9c6;" align="center">
The pipeline <strong>'.$pipeline.'</strong> has completed</p>';
	}
	

	$html_content .= '

<hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />

<table class="run-stats" style="border-spacing: 0; border-collapse: collapse; vertical-align: top; text-align: left; padding: 0;">
	<tr style="vertical-align: top; text-align: left; padding: 0;" align="left">
		<th style="text-align: right; padding-right: 10px; min-width: 70px;" align="right">
			Started
		</th>
		<td style="border-collapse: collapse !important; vertical-align: top; text-align: left; color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; line-height: 19px; font-size: 14px; margin: 0; padding: 0px 0px 10px;" align="left" valign="top">
			'.$startdate.'
		</td>
	</tr>
	<tr>
		<th style="text-align: right; padding-right: 10px; min-width: 70px;" align="right">
			Finished
		</th>
		<td style="border-collapse: collapse !important; vertical-align: top; text-align: left; color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; line-height: 19px; font-size: 14px; margin: 0; padding: 0px 0px 10px;" align="left" valign="top">
			'.$date.'
		</td>
	</tr>
	<tr>
		<th style="text-align: right; padding-right: 10px; min-width: 70px;" align="right">
			Duration
		</th>
		<td style="border-collapse: collapse !important; vertical-align: top; text-align: left; color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; line-height: 19px; font-size: 14px; margin: 0; padding: 0px 0px 10px;" align="left" valign="top">
			'.$duration.'
		</td>
	</tr>
	<tr>
		<th style="text-align: right; padding-right: 10px; min-width: 70px;" align="right">
			Working directory
		</th>
		<td style="border-collapse: collapse !important; vertical-align: top; text-align: left; color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; line-height: 19px; font-size: 14px; margin: 0; padding: 0px 0px 10px;" align="left" valign="top">
			<code style="font-family: \'Lucida Console\', Monaco, monospace; font-size: 12px; background: #efefef; padding: 3px 5px;">
				'.$cwd.'
			</code>
		</td>
	</tr>
</table>

<hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />

';

# Highlight any warnings
if($warnings > 0){
	$html_content .= '<h3 style="color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: left; line-height: 1.3; word-break: normal; font-size: 32px; margin: 0; padding: 0;" align="left">Custom Warnings Found</h3>';
	$html_content .= '<ul style="padding-left:20px;">';
	foreach my $file (sort keys %warninglines) {
		if(scalar @{$warninglines{$file}} > 0){
			$html_content .= '<li style="margin-bottom: 20px;"><span class="run-name" style="font-weight: bold; margin-bottom: 10px;">'.$file.'</span>';
			$html_content .= '<ul style="padding-left:20px;">';
			foreach my $warning (@{$warninglines{$file}}){
				$html_content .= '<li style="margin-top: 5px; background: #fcf8e3; color: #4D3E25; font-weight:bold; padding: 3px 5px;">'.$warning .'</li>';
			}
			$html_content .= '</ul></li>';
		}
	}
	$html_content .= '</ul><hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />';
}

# Highlight any highlights
if($highlights > 0){
	$html_content .= '<h3 style="color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: left; line-height: 1.3; word-break: normal; font-size: 32px; margin: 0; padding: 0;" align="left">Custom Highlights Found</h3>';
	$html_content .= '<ul style="padding-left:20px;">';
	foreach my $file (sort keys %highlightlines) {
		if(scalar @{$highlightlines{$file}} > 0){
			$html_content .= '<li style="margin-bottom: 20px;"><span class="run-name" style="font-weight: bold; margin-bottom: 10px;">'.$file.'</span>';
			$html_content .= '<ul style="padding-left:20px;">';
			foreach my $highlight (@{$highlightlines{$file}}){
				$html_content .= '<li style="margin-top: 5px;">'.$highlight .'</li>';
			}
			$html_content .= '</ul></li>';
		}
	}
	$html_content .= '</ul><hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />';
}

$html_content .= '

<h3 style="color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: left; line-height: 1.3; word-break: normal; font-size: 32px; margin: 0; padding: 0;" align="left">CF Status Messages</h3>
<ul class="status-messages" style="padding-left:20px;">';

foreach my $file (sort keys %cf_highlights) {
	$html_content .= '<li style="margin-bottom: 20px;">
	<span class="run-name" style="font-weight: bold; margin-bottom: 10px;">
	'.$file.'</span>';
	if(scalar @{$cf_highlights{$file}} > 0){
		$html_content .= '<ul style="padding-left:20px;">';
		foreach my $highlight (@{$cf_highlights{$file}}){
			$html_content .= '<li style="margin-top: 5px;';
			if($highlight =~ /error/i){
				$html_content .= ' background: #ebccd1; color: #222222; font-weight:bold; padding: 3px 5px;';
			}
			$html_content .= '">'.$highlight.'</li>';
		}
		$html_content .= '</ul>';
	}
	$html_content .= '</li>';
}

if(scalar @summarylines > 0){
	$html_content .= '<li style="margin-bottom: 20px;">
	<span class="run-name" style="font-weight: bold; margin-bottom: 10px;">
	Summary Modules</span><ul style="padding-left:20px;">';
    foreach my $ln (@summarylines){
        $html_content .= '<li style="margin-top: 5px;">'.$ln.'</li>';
    }
    $html_content .= '</ul></li>';
}

$html_content .= '
</ul>

<hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />
<h3 style="color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: left; line-height: 1.3; word-break: normal; font-size: 32px; margin: 0; padding: 0;" align="left">Commands Run</h3>

<ul class="status-messages commands-run" style="padding:0; list-style-type:none;">';

foreach my $file (sort keys %commands) {
	$html_content .= '<li style="margin-bottom: 20px;">
	<span class="run-name" style="font-weight: bold; margin-bottom: 10px;">
	'.$file.'</span>';
	if(scalar @{$commands{$file}} > 0){
		$html_content .= '<ul style="padding:0; list-style-type:none;">';
		foreach my $command (@{$commands{$file}}){
			$html_content .= '<li style="margin-top: 5px; font-family: \'Lucida Console\', Monaco, monospace; font-size: 12px; background: #efefef; padding: 3px 5px;">'.$command .'</li>';
		}
		$html_content .= '</ul>';
	}
	$html_content .= '</li>';
}

$html_content .= '</ul><hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />';
	
	
	
	
	
	
	#### SEND THE EMAIL
	my $to = $config{email};
	my $subject = "$pipeline pipeline complete";
	my $title = "Run Complete";
	
	if(CF::Helpers::send_email($subject, $to, $title, $html_content, $plain_content)){
		warn "Sent a pipeline e-mail notification to $to\n";
	} else {
		warn "Error! Problem whilst trying to send a pipeline e-mail notification to $to\n";
	}

	#### SAVE THE REPORTS
  my ($html_email, $text_email) = CF::Helpers::build_emails($title, $html_content, $plain_content);
	open (HTML,'>',$pipeline_id.'.html') or die "Can't write to $pipeline_id.html: $!";
	print HTML $html_email;
	close(HTML);
	open (PLAIN,'>',$pipeline_id.'.txt') or die "Can't write to $pipeline_id.txt: $!";
	print PLAIN $text_email;
	close(PLAIN);

} elsif($config{notifications}{complete} && !defined($config{email})){
	warn "Error! Tried to send run e-mail notification but no e-mail address found in config\n";
}


