#!/usr/bin/env python
"""
-----------------------
Preseq Plotting Module
-----------------------
Takes results from preseq and plots nice colourful complexity curves
using the plot_complexity_curves() function in the ngi_visualizations
Python package.

This CF module can be run either as a regular module: #plot_preseq
or as a pipeline summary module: >plot_preseq
If run as a summary module, all preseq curves will be plotted
in a single graph.

The ngi_visualizations package is available here:
https://github.com/ewels/ngi_visualizations

To install, run:
pip install git+https://github.com/NationalGenomicsInfrastructure/ngi_visualizations.git

 ******************************************************************
NB: Loading some modules (eg. cutadapt) can change your python path.
You can try running 'module load cutadapt' before the above pip installation.
If having permission errors, try 'pip install --user' instead.
 ******************************************************************

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
import sys
from CF import Helpers

# Module requirements
def mod_time_estimate(cf):
    num_files = cf['num_starting_merged_aligned_files']
    num_files = 1 if num_files < 1 else num_files
    return Helpers.minutes_to_timestamp(num_files * 6 * 60);

requirements = {
	'cores':   '1',
	'memory':  '3G',
	'modules': '',
	'time': mod_time_estimate
}

def make_preseq_plots(cf):
    """
    Main function to make plots using the plot_complexity_curves()
    function in the ngi_visualizations Python package.
    """
    # Load packages
    try:
        from ngi_visualizations.preseq_complexity_curves import plot_complexity_curves
    except ImportError, e:
        print("###CF Error: ngi_visualizations Python Package not installed.\n", file=sys.stderr)
        raise ImportError(e)

    # Are we running as a summary module?
    if 'summary_module' in cf['parameters']:
        timestart = datetime.datetime.now()

        # scrape the last file names from each run file
        preseq_files = []
        for cf in cf['runfns']:
            files, config = Helpers.parse_cf(cf, False, True)
            preseq_files.extend(files)

        # Make one single plot with all files
        plot_complexity_curves.plot_complexity_curves(preseq_files, output_name = "{}_preseq".format(cf['pipeline_name']))
        duration = str(datetime.datetime.now() - timestart)
        print("\n###CFCMD Ran summary ngi_visualizations.plot_complexity_curves.plot_complexity_curves()\n\n", file=sys.stderr)
        print("###CF Preseq summary plot successfully exited, took {}..\n".format(duration), file=sys.stderr)


    # Make one plot per file
    else:
        for fn in cf['files']:
            timestart = datetime.datetime.now()
            plot_complexity_curves.plot_complexity_curves(fn)
            duration = str(datetime.datetime.now() - timestart)
            print("\n###CFCMD Ran ngi_visualizations.plot_complexity_curves.plot_complexity_curves({})\n\n".format(fn), file=sys.stderr)
            print("###CF Preseq plot for {} successfully exited, took {}..\n".format(fn, duration), file=sys.stderr)


# Setup
if __name__ == "__main__":
    cf = Helpers.module_start(requirements, __doc__)
    make_preseq_plots(cf)
