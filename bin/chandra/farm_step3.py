#!/usr/bin/env python

"""
Check each candidate for hot pixels and nearest variable sources
"""

import argparse
import os
import sys

from XtDac.ChandraUtils import find_files
from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils.run_command import CommandRunner
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename
from XtDac.ChandraUtils.work_within_directory import work_within_directory

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Check for hot pixels/variable sources, add them to a '
                                                 'transient masterlist, and generate lightcurves')

    parser.add_argument("-o", "--obsid", help="Observation ID Numbers", type=int, required=True, nargs="+")
    parser.add_argument("-d", "--data_path", help="Path to directory containing data of all obsids", required=True,
                        type=str)
    parser.add_argument("-f", "--outfile", help="Name of output file which will contain filtered list of transients",
                        required=True)
    parser.add_argument("-s", "--resampleFactor",
                        help="Oversample the input image by this factor before processing",
                        type=int, default=5, required=False)

    # Get the logger
    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Get the command runner
    runner = CommandRunner(logger)

    args = parser.parse_args()

    data_path = sanitize_filename(args.data_path)

    with work_within_directory(data_path):

        for this_obsid in args.obsid:

            if not os.path.exists(str(this_obsid)):
                raise IOError("Directory not found for obsid %s" % this_obsid)

            with work_within_directory(str(this_obsid)):

                ccd_files = find_files.find_files('.', 'ccd*%s*_filtered.fits' % this_obsid)
                ccd_files = sorted(ccd_files)

                ccd_bb_files = find_files.find_files('.', 'ccd*%s*res.txt' % this_obsid)
                ccd_bb_files = sorted(ccd_bb_files)

                evtfile = find_files.find_files('.', '*%s*evt3.fits' % this_obsid)[0]

                if len(ccd_bb_files) != len(ccd_files):
                    raise RuntimeError("\n\nUnequal number of CCD files and BB files")

                for i in xrange(len(ccd_bb_files)):
                    og_file = os.path.basename(ccd_bb_files[i])

                    ccd_bb_file = ccd_bb_files[i]
                    print ccd_bb_file

                    ccd_file = ccd_files[i]
                    print ccd_file

                    check_hp_file = "check_hp_%s" % og_file

                    cmd_line = "xtc_remove_hot_pixels.py --obsid %s --evtfile %s --bbfile %s --outfile %s --debug no" \
                               % (this_obsid, ccd_file, ccd_bb_file, check_hp_file)

                    runner.run(cmd_line)

                    check_var_file = "check_var_%s" % og_file

                    cmd_line = "xtc_flag_near_variable_sources.py --bbfile %s --outfile %s --eventfile %s" \
                               % (check_hp_file, check_var_file, evtfile)

                    runner.run(cmd_line)

                check_var_files = find_files.find_files('.', 'check_var*%s*txt' % this_obsid)

            for check_var_file in check_var_files:
                cmd_line = "xtc_add_to_masterlist.py --bbfile %s --masterfile %s --evtfile %s" \
                           % (check_var_file, args.outfile, evtfile)

                runner.run(cmd_line)
