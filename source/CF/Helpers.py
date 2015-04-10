#!/usr/bin/env python
"""
Helpers.py
Module containing Cluster Flow helper functions written in Python
"""

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

from __future__ import print_function

import datetime
import os
import re
import sys
from inspect import stack


# Open the run file and parse it's contents.
# Input arguments:
#  - Run file path (required)
#  - Previous job id (optional)
# Returns a tuple containing:
#  - A list of file names
#  - A dictionary with config variables.
def parse_runfile (runfile_fn, prev_job_id=False, last_job=False):

    # Set up variables
    files = []
    config = {}
    config['notifications'] = {}
    config['references'] = {}
    comment_block = False
    last_jid = False

    # Read the file
    try:
        with open(runfile_fn) as runfile:
            for line in runfile:
                line = line.strip()

                # Ignore comment blocks
                if line[:2] == '/*':
                    comment_block = True
                    continue
                if line[:2] == '*/':
                    comment_block = False
                    continue

                # Get config variables
                if line[:1] == '@' and comment_block is False:
                    sections = line.split(None, 2)
                    cname = sections[0][1:]
                    if cname == 'notification':
                        config['notifications'][sections[1]] = 1
                    elif cname == 'reference':
                        ref_sections = sections[1].split(None, 2)
                        config['references'][ref_sections[0]] = ref_sections[1]
                    else:
                        config[cname] = sections[1]

                # Get files
                elif comment_block is False:
                    sections = line.split(None, 2)
                    try:
                        if last_job:
                            if sections[0] != last_jid:
                                prev_job_id = sections[0]
                                files = []
                            last_jid = sections[0]

                        if sections[0] == prev_job_id:
                            sections[1] = sections[1].strip()
                            files.append(sections[1])
                    except IndexError:
                        pass

    except IOError as e:
        print("Can't read run file: {}".format(runfile_fn))
        raise IOError(e)

    return (files, config)




# Used by modules when called. Either prints the help text, the module
# requirements or returns a hash with information for the run.
# Takes three variables - the @ARGV reference, a hash with keys for
# cores, memory, modules and time requirements (can be arrays, strings or
# subroutines) and the help text string.
def parse_runfile (reqs=False, helptext=False):

    # Find out where we're being called from
    script_fn = stack()[1][1]
    modname = script_fn.replace('.cfmod.py', '') # Always a python module
    modname = os.basename(modname)

    # Command line arguments
    parser = argparse.ArgumentParser("Cluster Flow Module {}".format(modname))
    parser.add_argument('--requirements', dest='get_requirements', action='store_true',
                        help="Request the cluster resources needed by the module")
    parser.add_argument('--run_fn', dest='runfns', type=str, action='append',
                        help="Path to the Cluster Flow run file(s) for this pipeline")
    parser.add_argument('--job_id', dest='job_id', type=str,
                        help="Cluster job ID for this job")
    parser.add_argument('--prev_job_id', dest='prev_job_id', type=str,
                        help="Cluster job ID for the previous job")
    parser.add_argument('--cores', dest='cores', type=int,
                        help="Number of cores to be used by the module")
    parser.add_argument('--mem', dest='mem', type=str,
                        help="Amount of memory to be used by the module")
    parser.add_argument('--param', dest='params_list', type=str, action='append',
                        help="Extra parameters to be used")
    parser.add_argument('--help', dest='show_help', action='store_true',
                        help="Print help")

    kwargs = vars(parser.parse_args())

    # Variable massaging
    params = {}
    for param in kwargs.params_list:
        try:
            (key, var) = param.split('=', 1)
            params[key] = var
        except ValueError:
            params[param] = True

    # Initialise the run file dict
    runfile = {
        'run_fn':       kwargs.runfns[0],
        'run_fns':      kwargs.runfns,
        'modname':      modname,
        'mod_fn':       script_fn,
        'job_id':       kwargs.job_id,
        'prev_job_id':  kwargs.prev_job_id,
        'cores':        kwargs.cores,
        'memory':       kwargs.mem,
        'params':       params
    }

    # Print help
    if show_help:
        if helptext:
            print helptext
        else:
            print "No help text for this module"
        quit()
    }

    # Check that we have a run file, need it to go any further
    assert len(kwargs.runfns[0]) > 0

    ###
    # Module requirements
    ###
    if kwargs.get_requirements:

        # Parse the run file(s)
        parse_runfile(runfile)

        # Run through the supplied requirements and print
        for key in reqs.keys():

            # Been given a function
            if type(reqs[key]) is 'function':
                print "{}: {}\n".function(key, reqs[key](runfile))

            # Been given an array
            if type(reqs[key]) is 'list':
                if key is 'cores':
                    print "cores: {}\n".format(allocate_cores(kwargs.cores, int(reqs[key][0]), int(reqs[key][1])))
                elif key is 'memory':
                    print "memory: {}\n".bytes_to_human_readable(allocate_memory(kwargs.mem, reqs[key][0], reqs[key][1]))
                else:
                    arrstring = ', '.join(reqs[key])
                    print "{}: {}\n".format(key, arrstring)


            # Been given a string
            elif(type(reqs[key]) is 'str'){
                print "{}: {}\n".format(key, reqs[key])

        # Printed everything - exit
        quit()



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

































# Used by modules at run time. Parses the incoming command line arguments
# and run file. Returns a dictionary with lots of information.
def load_runfile_params (params):

    # Pop the known variables off the front of the list
    try:
        runfile = params.pop(0)
        job_id = params.pop(0)
        prev_job_id = params.pop(0)
        cores = params.pop(0)
        mem = params.pop(0)
        # params should now contain any remaining paramters
    except IndexError as e:
        raise KeyError("Error: Required parameters not supplied to load_runfile_params() in Helpers.py")

    # Print the module header
    if 'hide_log_header' not in params:
        modname = __file__
        if modname[-6:] == '.cfmod':
            modname = modname[:-6]
        if len(modname) > 0:
            modname = "Module:\t\t\t{}\n".format(modname);

        paramterstring = ', '.join(params)
        date = datetime.datetime.now().strftime("%H:%M, %d-%m-%Y")
        dashes = '-' * 80
        print ("\n{dashes}\n{modname}Run File:\t\t{runfile}\nJob ID:\t\t\t{job_id}\nPrevious Job ID:\t{prev_job_id}\nParameters:\t\t{paramterstring}\nDate & Time:\t\t{date}\n{dashes}\n\n".format(**locals()), file=sys.stderr);

    files, config = parse_runfile(runfile, prev_job_id)

    # If we don't have any input files, bail now
    if len(files) == 0 and prev_job_id != 'null':
        # Allow this for summary modules
        if 'summary_module' not in params:
            print("\n###CF Error! No file names found from job {}. Exiting...\n\n".format(prev_job_id, file=sys.stderr))
            raise IOError

    returndict = {
        'files': files,
        'runfile': runfile,
        'job_id': job_id,
        'prev_job_id': prev_job_id,
        'cores': cores,
        'mem': mem,
        'parameters': params,
        'config': config
    }

    return returndict
