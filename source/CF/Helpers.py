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

import argparse
import datetime
import inspect
import math
import os
import re
import sys


# Used by modules when called. Either prints the help text, the module
# requirements or returns a hash with information for the run.
# Takes three variables - the @ARGV reference, a hash with keys for
# cores, memory, modules and time requirements (can be arrays, strings or
# subroutines) and the help text string.
def module_start (reqs=False, helptext=False):

    # Find out where we're being called from
    script_fn = inspect.stack()[1][1]
    modname = script_fn.replace('.cfmod.py', '') # Always a python module
    modname = os.path.basename(modname)

    # Command line arguments
    parser = argparse.ArgumentParser("Cluster Flow Module {}\n\n{}".format(modname, helptext))
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

    kwargs = vars(parser.parse_args())

    # Variable massaging
    params = {}
    if kwargs['params_list'] is not None:
        for param in kwargs['params_list']:
            try:
                (key, var) = param.split('=', 1)
                params[key] = var
            except ValueError:
                params[param] = True

    # Check that we have at least run file, need one to go any further
    assert kwargs['runfns'] is not None and len(kwargs['runfns']) > 0

    # Initialise the run file dict
    cf = {
        'run_fn':       kwargs['runfns'][0],
        'run_fns':      kwargs['runfns'],
        'modname':      modname,
        'mod_fn':       script_fn,
        'job_id':       kwargs['job_id'],
        'prev_job_id':  kwargs['prev_job_id'],
        'cores':        kwargs['cores'],
        'memory':       kwargs['mem'],
        'params':       params
    }

    # Parse things from the job ID
    if kwargs['job_id']:
        jobid_re = re.match('^cf_(.+)_(\d{10})_(.+)_\d{1,3}$', kwargs['job_id'])
        if jobid_re:
            cf['pipeline_name'] = jobid_re.group(1)
            cf['pipeline_started'] = jobid_re.group(2)
            cf['pipeline_id'] = "{}_{}".format(cf['pipeline_name'], cf['pipeline_started'])

    ###
    # Module requirements
    ###
    if kwargs['get_requirements']:

        # Parse the run file(s)
        parse_runfile(cf)

        # Run through the supplied requirements and print
        for key in reqs.keys():

            # Been given a function
            if hasattr(reqs[key], '__call__'):
                print ("{}: {}".format(key, reqs[key](cf)))

            # Been given an array
            if type(reqs[key]) is list:
                if key == 'cores':
                    print ("cores: {}".format(allocate_cores(kwargs['cores'], int(reqs[key][0]), int(reqs[key][1]))), file=sys.stout)
                elif key == 'memory':
                    print ("memory: {}".bytes_to_human_readable(allocate_memory(kwargs['mem'], reqs[key][0], reqs[key][1])), file=sys.stout)
                else:
                    print ("{}: {}".format(key, ', '.join(reqs[key])), file=sys.stout)


            # Been given a string
            elif type(reqs[key]) is str:
                print ("{}: {}".format(key, reqs[key]))

        # Printed everything - exit
        quit()



    ###
    # Running for real
    ###
    # check that we have extra parameters that we need
    if not kwargs['job_id'] or not kwargs['prev_job_id']:
        print ("###CF Error: Need job ID and previous job ID! Missing for module {}\n".format(modname))
        raise TypeError

    if not kwargs['cores'] or not kwargs['mem']:
        print ("###CF Error: Need allocated cores and memory! Missing for module {}\n".format(modname))
        raise TypeError


    # Parse everything we can from the run file
    parse_runfile(cf)

    # Print the module header if needed
    if 'hide_log_header' in params:

        dashes = '-' * 80
        header = "\n{}\nModule:\t\t\t{}\n".format(dashes, modname)

        if 'summary_module' in params:
            header += "Summary Module:\t\tYes\n"

        header += "Run File:\t\t{}\nJob ID:\t\t\t{}\nPrevious Job ID:\t{}\n".format(cf['run_fn'], cf['job_id'], cf['prev_job_id'])

        if len(params) > 0:
            s_params = ['{} = {}'.format(k,v) for k,v in params.iteritems()]
            header += "Parameters:\t\t{}\n".format(', '.join(s_params))

        datestring = datetime.datetime.now().strftime("%H:%M, %d-%m-%Y")
        header += "Date & Time:\t\t{}\n{}\n\n".format(datestring, dashes)

        print (header, file=sys.stderr)


    # If we don't have any input files, bail now
    if ('prev_job_files' not in cf or len(cf['prev_job_files']) == 0) \
        and (kwargs['prev_job_id'] and kwargs['prev_job_id'] != 'null' and 'summary_module' not in params):
            print ("\n###CF Error! No file names found from job {}. Exiting...\n\n".format(kwargs['prev_job_id']))
            raise TypeError

    # Return the dictionary of dreams
    return cf



# Open the run file and parse it's contents.
# Input arguments: cf dict
# Returns: Nothing (updates dictionary in place)
def parse_runfile (cf):

    # Set up new dict variables if we need them
    if 'refs' not in cf:
        cf['refs'] = {}
    if 'config' not in cf:
        cf['config'] = {}
    if 'notifications' not in cf['config']:
        cf['config']['notifications'] = {}
    if 'prev_job_files' not in cf:
        cf['prev_job_files'] = []
    if 'starting_files' not in cf:
        cf['starting_files'] = []
    if 'files' not in cf:
        cf['files'] = {}
    if 'num_starting_files' not in cf:
        cf['num_starting_files'] = 0
    if 'num_starting_merged_files' not in cf:
        cf['num_starting_merged_files'] = 0
    if 'num_starting_merged_aligned_files' not in cf:
        cf['num_starting_merged_aligned_files'] = 0

    # Go through each run file
    for runfn in cf['run_fns']:
        comment_block = False
        try:
            with open(runfn) as runfile:
                for line in runfile:

                    line = line.strip()

                    # Ignore comment blocks
                    if line[:2] == '/*':
                        comment_block = True
                        continue
                    if line[:2] == '*/':
                        comment_block = False
                        continue

                    # Helper stuff
                    if line[:9] == 'Pipeline:':
                        cf['pipeline_name'] = line[10:]

                    # Get config variables
                    elif line[:1] == '@' and comment_block is False:
                        sections = line.split(None, 2)
                        cname = sections[0][1:]
                        if not sections[1]:
                            sections[1] = 1
                        if cname == 'notification':
                            cf['config']['notifications'][sections[1]] = 1
                        elif cname == 'reference':
                            ref_sections = sections[1].split(None, 2)
                            cf['refs'][ref_sections[0]] = ref_sections[1]
                        else:
                            cf['config'][cname] = sections[1]

                    # Note - could get the pipeline tree here. Currently no need.

                    # Get files
                    elif not comment_block is False and line[:1] not in ['@', '#', '>']:
                        sections = line.split(None, 2)

                        # Previous job files
                        if 'prev_job_id' in cf and sections[0] == cf['prev_job_id']:
                            cf['prev_job_files'].append(sections[1])

                        # Starting files
                        if sections[0] == 'start_000':
                            cf['starting_files'].append(sections[1])
                            cf['num_starting_files'] += 1

                        # All files, by module
                        if sections[0] not in cf['files']:
                            cf['files'][sections[0]] = [sections[1]]
                        else:
                            cf['files'][sections[0]].append(sections[1])

        except IOError as e:
            print("Can't read run file: {}".format(runfn))
            raise IOError(e)

    # Figure out how many merged files we're likely to have
    regex = False
    if 'merge_regex' in cf['config'] and len(cf['config']['merge_regex']) > 0:
        regex = cf['config']['merge_regex']
    if 'regex' in cf['params']:
        regex = cf['params']['regex']

    if regex:
        file_sets = {}
        for file in cf['starting_files']:
            m = re.search(regex, file)
            try:
                file_sets[m.groups()[0]] = 1
            except AttributeError:
                file_sets[file] = 1

        cf['num_starting_merged_files'] = len(file_sets.keys())
    else:
        cf['num_starting_merged_files'] = cf['num_starting_files']

    # How many aligned files? Just check for paired end or single end config
    if 'force_paired_end' in cf['config']:
        cf['num_starting_merged_aligned_files'] = cf['num_starting_merged_files'] / 2
    else:
        cf['num_starting_merged_aligned_files'] = cf['num_starting_merged_files']

    # No need to return anything, as dict updates in place




def timestamp_to_minutes(timestamp):
    """
    Function to convert a SLURM style timestamp to minutes.
    Acceptable input formats:
        minutes
        minutes:seconds
        hours:minutes:seconds
        days-hours
        days-hours:minutes
        days-hours:minutes:seconds
    """
    days = 0
    hours = 0
    minutes = 0
    # seconds = 0
    re_hh = "([0-2]?[0-9])"
    re_mm = "([0-5]?[0-9])"

    # Make sure that we don't have any illegal characters
    if re.match("[^\d\-\:]", timestamp):
        return 0


    # First - days
    days_re = re.match('^(\d+)\-', timestamp)
    if days_re:
        days = days_re.group(1)
        timestamp = re.sub('^(\d+)\-', '', timestamp)

    # Progressively match colons to figure out what the rest is
    hms_re = re.match('^{}\:{}\:{}$'.format(re_hh, re_mm, re_mm), timestamp)
    if hms_re:
        hours = hms_re.group(1)
        minutes = hms_re.group(2)
        # seconds = hms_re.group(3)
    else:
        if days > 0:
            hm_re = re.match('^{}\:{}$'.format(re_hh, re_mm), timestamp)
            if hm_re:
                hours = hm_re.group(1)
                minutes = hm_re.group(2)
            else:
                h_re = re.match('^{}$'.format(re_hh), timestamp)
                if h_re:
                    hours = h_re.group(1)
        else:
            ms_re = re.match('^{}\:{}$'.format(re_mm, re_mm), timestamp)
            if ms_re:
                minutes = ms_re.group(1)
                # seconds = ms_re.group(2)
            else:
                m_re = re.match('^{}$'.format(re_hh), timestamp)
                if m_re:
                    minutes = m_re.group(1)

    total_minutes = 0
    total_minutes += days * 24 * 60
    total_minutes += hours * 60
    total_minutes += minutes
    # We ignore seconds. Seriously, who does that?

    return total_minutes


# Function to convert minutes to a SLURM time stamp. Reverse of above.
def minutes_to_timestamp(minutes):

    minutes = int(minutes)

    days = int(minutes / (24 * 60))
    minutes -= days * 24 * 60
    minutes = 0 if minutes < 0 else minutes

    hours = int(minutes / 60)
    minutes -= hours * 60
    minutes = 0 if minutes < 0 else minutes

    minutes = "{0:02d}".format(minutes)

    if days > 0:
        return "{}-{}:{}".format(days, hours, minutes)
    elif hours > 0:
        return "{}:{}:00".format(hours, minutes)
    else:
        return minutes



# Function to take human readable memory string and return bytes
def human_readable_to_bytes(memory):

    # Find the suffix if we have one
    suffix = ''
    match = re.match('([tgmk])b?', memory.lower())
    if match:
        suffix = match.group(1)

    # Remove non-numeric (except decimal places)
	memory = re.sub('[^\d\.]', '', memory)

	if suffix == 't':
		return memory * 1000000000000
	elif suffix == 'g':
		return memory * 1000000000
	elif suffix == 'm':
		return memory * 1000000
	elif suffix == 'k':
		return memory * 1000


# Function to take bytes and return a human readable memory string
def bytes_to_human_readable (mybytes):

	# Remove non-numeric
	mybytes = re.sub('\D', '', mybytes)

	if int(mybytes/1000000000) > 0:
		return "{}G".format(math.ceil(mybytes/1000000000))
	elif int(mybytes/1000000) > 0:
		return "{}M".format(math.ceil(mybytes/1000000))
	elif int(mybytes/1000) > 0:
		return "{}K".format(math.ceil(mybytes/1000))
	else:
		return "{}B".format(mybytes)

# Simple function to take bytes and return a human readable memory string
# NB: Rounds up with math.ceil()
def mem_return_mbs(mem):
    mem = human_readable_to_bytes(mem)
    return math.ceil(mem/1000000);


# Take allocated cores, minimum, maximum and return best value
def allocate_cores(allocated, min_cores, max_cores):

	if allocated > max_cores:
		return max_cores
	elif allocated < min_cores:
		return min_cores
	else:
		return allocated

# Take allocated memory, minimum, maximum and return best value
def allocate_memory(allocated, min_mem, max_mem):

    allocated = human_readable_to_bytes(allocated);
    max_mem = human_readable_to_bytes(max_mem);
    min_mem = human_readable_to_bytes(min_mem);

    if allocated > max_mem:
        return max_mem
    elif allocated < min_mem:
        return min_mem
    else:
        return allocated
