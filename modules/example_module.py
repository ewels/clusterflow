#!/usr/bin/env python
"""
example_module.py

This is an example module for Cluster Flow.
It does not have any real functionality and should not be used in
production. To avoid its accidental use in Cluster Flow, the 
file name does not end in .cfmod
If developing a new module from this example, remember to change
the filename to MODULENAME.cfmod(.py)

Note: This code is non-functional, so difficult to test. As such,
it's basically been written blind. Apologies if there are any really
horrible errors.
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
from CF import Helpers

import argparse
import datetime
import os
import shlex
import subprocess
import sys

from my_favourite_python_package import imported_functions

def do_main_module_function(parameters, required_cores=False, required_mem=False, required_modules=False, runfn=None, print_help=False):
    """
    -----------------------
    Example Python Module
    -----------------------
    Takes results from preseq and plots nice colourful complexity curves
    using the plot_complexity_curves() function in the ngi_visualizations
    Python package. This package is available here:
    https://github.com/ewels/ngi_visualizations
    """
    
    #
    # JOB INITIALISATION
    # These flags are used on the head node to allocate resources
    #

    # --cores specifies how many cores are recommended
    # Return how many cores we need. This can be more or less than
    # the number recommended.
    if required_cores:
        print ('1', file=sys.stdout)
        sys.exit(0)

    # --mem specifies how much memory is recommended
    # Return how much we need. This can be more or less than the amount recommended.
    if required_mem:
        print ('3G', file=sys.stdout)
        sys.exit(0)

    # --modules. Return comma seperated names of any
    # environment modules which should be loaded.
    if required_modules:
        print ('example_module,samtools', file=sys.stdout)
        sys.exit(0)

    # --help. Print help. This code prints the function docstring above
    if print_help:
        print (make_preseq_plots.__doc__, file=sys.stdout)
        sys.exit(0)
    

    #
    # MODULE CODE
    # If we get this far, the module is being run in the pipeline
    #

    # Start the clock
    timestart = datetime.datetime.now()

    # Parse the incoming parameters and the run file
    # This function returns a number of configuration variables as well 
    # as the incoming filenames as a dictionary with the following keys:
    #  files
    #  runfile
    #  job_id
    #  prev_job_id
    #  cores
    #  mem
    #  parameters
    #  config
    p = Helpers.load_runfile_params(parameters)
    
    #
    # RUNNING PYTHON FUNCTIONS
    # This block of code serves as an example of how to execute imported Python code
    #

    # Run the imported python function with the list of filenames
    output_files = imported_functions.my_def(p['files'])
    print("\n###CFCMD Ran my_def() from my_favourite_python_package.imported_functions\n\n", file=sys.stderr)

    # How long did it take?
    duration = str(datetime.datetime.now() - timestart)
    print("###CF Example python command successfully exited, took {}..\n".format(duration), file=sys.stderr)

    # Write the output filen ames to the run file
    for output_fn in output_files:
        # Check we can find our output file
        if os.path.isfile(output_fn):
            # Print the current job ID and the output filename to the run file
            # This is so that subsequent modules can use this output
            try:
                with open(p['runfile'], 'a') as runfile:
                    print("{}\t{}\n".format(p['job_id'], output_fn), file=runfile)
            except IOError as e:
                print("###CF Error: Can't write to {}\n".format(runfile))
                raise IOError(e)



    #
    # RUNNING EXTERNAL SYSTEM COMMANDS
    # This block of code is an example of how to call external programs
    # Obviously, this will over-write the output_files variable from above
    #

    # Print version information about the program to be executed if we can
    print("---------- < module > version information ----------\n", file=sys.stderr)
    print(subprocess.check_output(shlex.split('MY_COMMAND --version'), file=sys.stderr))  
    print("\n------- End of < module > version information ------\n", file=sys.stderr)
    
    # NOTE - missing
    # This code should really take the list of input files and split them into
    # paired end and single end files, but I haven't written the helper functions
    # in Python to do this yet. Raise an issue on GitHub if you need this.

    # Loop through the files
    for fn in p['files']:

        # What's our output filename?
        output_fn = "{}_processed.output".format(fn)

        # Put the command together
        cmd = "my_command -c {} -m {} -g {} -i {} -o {}".format(p['cores'], p['mem'], p['references']['fasta'], fn, output_fn)
        print("\n###CFCMD {}\n\n".format(cmd), file=sys.stderr)

        # Run the system command
        if subprocess.call(shlex.split(cmd)) == 0:

            # How long did it take?
            duration = str(datetime.datetime.now() - timestart)
            
            # Print a success message to the log file which will be e-mailed out
            print("###CF Example module system call was successful,  {}..\n".format(duration), file=sys.stderr)

            # Check we can find our output file
            if os.path.isfile(output_fn):

                # Print the current job ID and the output filename to the run file
                # This is so that subsequent modules can use this output
                try:
                    with open(p['runfile'], 'a') as runfile:
                        print("{}\t{}\n".format(p['job_id'], output_fn), file=runfile)
                except IOError as e:
                    print("###CF Error: Can't write to {}\n".format(runfile))
                    raise IOError(e)

            # Oops - can't find the output file!
            else:
                print("###CF Error: Example module output file {} not found\n".format(output_fn), file=sys.stderr)
                raise IOError

        # Command returned a non-zero code, something went wrong
        else:
            print("###CF Error: Example module failed for input file {}\n".format(fn), file=sys.stderr)
            raise SystemError
    
    
    

if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("Make a scatter plot of FPKM counts between conditions")
    parser.add_argument('--cores', dest='required_cores', action='store_true',
                        help="Request the number of cores needed by the module.")
    parser.add_argument('--mem', dest='required_mem', action='store_true',
                        help="Request the amount of memory needed by the module.")
    parser.add_argument('--modules', dest='required_modules', action='store_true',
                        help="Request the names of environment modules needed by the module.")
    parser.add_argument('--runfn', dest='runfn', type=str, default=None,
                        help="Path to the Cluster Flow run file for this pipeline")
    parser.add_argument('parameters', nargs='*', help="List of parameters.")
    kwargs = vars(parser.parse_args())
    
    # Call do_main_module_function()
    do_main_module_function(**kwargs)


