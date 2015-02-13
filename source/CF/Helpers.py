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
import re
import sys


# Open the run file and parse it's contents.
# Input arguments:
#  - Run file path (required)
#  - Previous job id (optional)
# Returns a tuple containing:
#  - A list of file names
#  - A dictionary with config variables.
def parse_runfile (runfile_fn, prev_job_id=False):
    
    # Set up variables
    files = []
    config = {}
    config['notifications'] = {}
    config['references'] = {}
    comment_block = False
    
    # Read the file
    try:
        with open(runfile_fn) as runfile:
            for line in fh:
                line = line.strip()
                
                # Ignore comment blocks
                if line[:2] == '/*':
                    comment_block = True
                    continue
                if line[:2] == '*\\':
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
                    if sections[0] == prev_job_id:
                        sections[1] = sections[1].strip()
                        files.append(sections[1])
    
    except IOError as e:
        print("Can't read run file: {}".format(runfile_fn))
        raise IOError(e)

    return (files, config)







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
        print ("\n{dashes}\n{modname}Run File:\t\t{runfile}\nJob ID:\t\t\t{job_id}\nPrevious Job ID:\t{prev_job_id}\nParameters:\t\t{paramterstring}\nDate & Time:\t\t{date}\n{dashes}\n\n".format(dashes, modname, runfile, job_id, prev_job_id, paramterstring, date), file=sys.stderr);
    
    files, config = parse_runfile(runfile, prev_job_id)
    
    # If we don't have any input files, bail now
    if len(files) == 0 and prev_job_id != 'null':
        print("\n###CF Error! No file names found from job {}. Exiting...\n\n".format(prev_job_id, file=sys.stderr))
        raise IOError
    
    returndict = {
        'files': files,
        'runfile': runfile,
        'job_id': job_id,
        'prev_job_id': prev_job_id,
        'cores': cores,
        'mem': mem,
        'parameters': parameters,
        'config': config
    }
    return returndict