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

# Module requirements
my %requirements = (
	'cores' 	=> '1',
	'memory' 	=> '1G',
	'modules' 	=> '',
	'time' 		=> '10'
);

# Help text
my $helptext = "\nThis is a core module which is executed when a single run has finished.\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);


# MODULE
my $pipeline = $cf{'pipeline_name'};
my $outfn = $cf{'params'}{'outfn'};

# Find current directory
my $cwd = getcwd()."/";

# Parse the outfile (the run logs)
my @LOG_HIGHLIGHT_STRINGS = @CF::Constants::LOG_HIGHLIGHT_STRINGS;
my @LOG_WARNING_STRINGS = @CF::Constants::LOG_WARNING_STRINGS;
my @cf_highlights;
my @commands;
my @warninglines;
my @highlightlines;
my %modules;
my %modkeys;
my $unrecognised = "";
my $i = 1;
my $errors = 0;
my $warnings = 0;
my $highlights = 0;
open (IN,'<',$outfn) or die "Can't read ".getcwd()."/$outfn - $!";
while(<IN>){

	chomp;

	# Ignore crap
	if(/^Warning: no access to tty/ || /^Thus no job control in this shell/ || /sh: module: line 1: syntax error: unexpected end of file/ || /sh: error importing function definition for/){
		next;
	}

	# First, strip the module identifier
	if(s/^###CF_(.*?)://){
		$modules{$1} .= "$_\n";
		$modkeys{$1} = $i;
	} else {
		$unrecognised .= "$_\n";
	}

	# Commands run
	if(/^###CFCMD/){
		push (@commands, substr($_, 9));

	# Highlight statuses
	} elsif(/^###CF/){
		push (@cf_highlights, substr($_, 6));
	} else {
		# Count any custom string highlights
		foreach my $highlight_string (@LOG_HIGHLIGHT_STRINGS){
			if(/$highlight_string/i){
				$highlights++;
				push (@highlightlines, $_);
			}
		}

		# Count any custom string errors
		foreach my $warning_string (@LOG_WARNING_STRINGS){
			if(/$warning_string/i){
				$warnings++;
				push (@warninglines, $_);
			}
		}
	}

	# Count out any CF errors
	if(/error/i){
		if (/^###CF/){
			$errors++;
		}
	}

	# Increment counter for module keys hash
	$i++;

}
close (IN);

# Make new, cleaned output
my $outfile = "";
# Go through module keys hash sorting by values (the last counter value seen)
foreach my $key (sort { $modkeys{$a} <=> $modkeys{$b} or $a cmp $b } keys %modkeys) {
	$outfile .= $modules{$key}."\n\n\n";
}
# count non-whitespace characters only
if($unrecognised =~ m/\S/g){
	$outfile .= ("-"x80)."\nUnrecognised output: (indicated with \"-> <-\")\n".("-"x80)."\n";
    $outfile .= "-> ".join(" <-\n-> ", split "\n", $unrecognised)." <-\n\n";
}
$outfile .= ("="x80)."\n\n\n";

# Print run finish status to outfile
my $date = strftime ("%H:%M %d-%m-%Y", localtime);
$outfile .= "\n\n\n###CF Run finished at $date\n\n";

# Write out the cleaned runfile without any crap
open (OUT,'>',$outfn) or die "Can't write to $outfn: $!";
print OUT $outfile;
close OUT;

# Send e-mail to submitter, if the config demands it
if($cf{'config'}{'notifications'}{'run'} && defined($cf{'config'}{'email'})){

	# Write the plain-text e-mail body
my $plain_content = "A run in the pipeline $pipeline has completed";
if($errors > 0){
	$plain_content .= " **with errors**";
} elsif ($warnings > 0){
	$plain_content .= " **with warnings**";
}
$plain_content .= ".

Finished: $date
Working Directory: $cwd

";

if($warnings > 0){
	$plain_content .= "\n\n\n\n\n
===========================
== Custom Warnings Found ==
===========================
- ".join("\n - ", @warninglines)."\n\n\n\n\n";
}
if($highlights > 0){
	$plain_content .= "\n\n\n\n\n
===========================
== Custom Highlights Found ==
===========================
- ".join("\n - ", @highlightlines)."\n\n\n\n\n";
}

$plain_content .= ".
========================
== CF status messages ==
========================
- ".join("\n- ", @cf_highlights)."
\n\n\n\n\n
==================
== Commands Run ==
==================
- ".join("\n", @commands)."

";

$plain_content .= "

== Full log output ==
".$outfile."\n\n\n";






	# Write the html e-mail body
	# Inline styles make me want to stab my eyes out, but we have to do this for
	# e-mail readers such as Gmail which strip any header CSS.
	my $html_content = '';
	if($errors > 0){
		$html_content .= '<p class="completion-leader" style="color: #723736; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: center; line-height: 19px; font-size: 18px; border-radius: 5px; background: #f2dede; margin: 0 0 15px; padding: 10px; border: 1px solid #ebccd1;" align="center">
A <strong>run</strong> in the pipeline <strong>'.$pipeline.'</strong> has completed <strong>with errors</strong></p>';
	} elsif($warnings > 0){
		$html_content .= '<p class="completion-leader" style="color: #4D3E25; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: center; line-height: 19px; font-size: 18px; border-radius: 5px; background: #fcf8e3; margin: 0 0 15px; padding: 10px; border: 1px solid #faebcc;" align="center">
The pipeline <strong>'.$pipeline.'</strong> has completed <strong>with warnings</strong></p>';
	} else {
		$html_content .= '<p class="completion-leader" style="color: #0A440B; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: center; line-height: 19px; font-size: 18px; border-radius: 5px; background: #dff0d8; margin: 0 0 15px; padding: 10px; border: 1px solid #d6e9c6;" align="center">
A <strong>run</strong> in the pipeline <strong>'.$pipeline.'</strong> has completed</p>';
	}
	$html_content .= '

<hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />

<table class="run-stats" style="border-spacing: 0; border-collapse: collapse; vertical-align: top; text-align: left; padding: 0;">
	<tr style="vertical-align: top; text-align: left; padding: 0;" align="left">
		<th style="text-align: right; padding-right: 10px; min-width: 70px;" align="right">
			Started
		</th>
		<td style="border-collapse: collapse !important; vertical-align: top; text-align: left; color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; line-height: 19px; font-size: 14px; margin: 0; padding: 0px 0px 10px;" align="left" valign="top">
			'.$date.'
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
	foreach my $warning (@warninglines){
		$html_content .= '<li style="margin-top: 5px; background: #fcf8e3; color: #4D3E25; font-weight:bold; padding: 3px 5px;">'.$warning.'</li>';
	}
	$html_content .= '</ul><hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />';
}

# Highlight any highlights
if($highlights > 0){
	$html_content .= '<h3 style="color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: left; line-height: 1.3; word-break: normal; font-size: 32px; margin: 0; padding: 0;" align="left">Custom Highlights Found</h3>';
	$html_content .= '<ul style="padding-left:20px;">';
	foreach my $highlight (@highlightlines){
		$html_content .= '<li style="margin-top: 5px;">'.$highlight.'</li>';
	}
	$html_content .= '</ul><hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />';
}

$html_content .= '
<h3 style="color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: left; line-height: 1.3; word-break: normal; font-size: 32px; margin: 0; padding: 0;" align="left">CF Status Messages</h3>
<ul class="status-messages" style="padding-left:20px;">';
$html_content .= '<ul style="padding-left:20px;">';
foreach my $highlight (@cf_highlights){
	$html_content .= '<li style="margin-top: 5px;';
	if($highlight =~ /error/i){
		$html_content .= ' background: #ebccd1; color: #222222; font-weight:bold; padding: 3px 5px;';
	}
	$html_content .= '">'.$highlight.'</li>';
}
$html_content .= '
</ul>

<hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />
<h3 style="color: #222222; font-family: \'Helvetica\', \'Arial\', sans-serif; font-weight: normal; text-align: left; line-height: 1.3; word-break: normal; font-size: 32px; margin: 0; padding: 0;" align="left">Commands Run</h3>

<ul class="status-messages commands-run" style="padding:0; list-style-type:none;">
	<li style="margin-top: 5px; font-family: \'Lucida Console\', Monaco, monospace; font-size: 12px; background: #efefef; padding: 3px 5px;">'.
	join('</li><li style="margin-top: 5px; font-family: \'Lucida Console\', Monaco, monospace; font-size: 12px; background: #efefef; padding: 3px 5px;">', @commands)
	.'</li>
</ul>

<hr style="color: #d9d9d9; height: 1px; background: #d9d9d9; border: none;" />';







	#### SEND THE EMAIL
	my $to = $cf{'config'}{'email'};
	my $subject = "$pipeline run compete - $cf{run_fn}";
	my $title = "Run Complete";

	if(CF::Helpers::send_email($subject, $to, $title, $html_content, $plain_content)){
		warn "Sent a pipeline e-mail notification to $to\n";
	} else {
		warn "Error! Problem whilst trying to send a pipeline e-mail notification to $to\n";
	}

} elsif($cf{'config'}{'notifications'}{'complete'} && !defined($cf{'config'}{'email'})){
	warn "Error! Tried to send run e-mail notification but no e-mail address found in config\n";
}
