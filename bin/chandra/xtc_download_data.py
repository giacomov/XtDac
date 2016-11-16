#!/usr/bin/env python

"""
Download files by observation ID
"""

import argparse
import os
import sys

from XtDac.ChandraUtils import find_files
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename
from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils.data_package import DataPackage
from XtDac.ChandraUtils.run_command import CommandRunner
from XtDac.ChandraUtils.work_within_directory import work_within_directory
from XtDac.ChandraUtils.configuration import get_configuration

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run following steps: download files by obsid')

    parser.add_argument("-c", "--config_file", help="Path to the configuration file", type=str, required=True)

    parser.add_argument("-o", "--obsid", help="Observation ID Numbers", type=int, required=True, nargs='+')

    args = parser.parse_args()

    # Get the logger
    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Get the command runner
    runner = CommandRunner(logger)

    # Get the configuration
    config = get_configuration(args.config_file)

    # Sanitize the workdir
    data_repository = sanitize_filename(config['data repository'])
    region_repository = sanitize_filename(config['region repository'])

    with work_within_directory(data_repository):

        # Download files

        for this_obsid in args.obsid:

            regdir_this_obsid = os.path.join(region_repository, str(this_obsid))

            if os.path.exists(regdir_this_obsid):

                cmd_line = "xtc_download_by_obsid.py --obsid %d" % this_obsid

                runner.run(cmd_line)

                try:

                    evtfile = os.path.basename(find_files.find_files(os.getcwd(), '*%s*evt3.fits' % this_obsid)[0])
                    tsvfile = os.path.basename(find_files.find_files(os.getcwd(), "%s.tsv" % this_obsid)[0])
                    expmap = os.path.basename(find_files.find_files(os.getcwd(), "*%s*exp3.fits*" % this_obsid)[0])
                    fov = os.path.basename(find_files.find_files(os.getcwd(), "*%s*fov3.fits*" % this_obsid)[0])

                except IndexError:

                    raise RuntimeError("\n\n\nCould not find one of the downloaded files for obsid %s. Exiting..."
                                       % this_obsid)

                data_package = DataPackage(str(this_obsid), create=True)

                data_package.store("evt3", evtfile, "Event file (Level 3) from the CSC", move=True)
                data_package.store("tsv", tsvfile, "TSV file from the CSC", move=True)
                data_package.store("exp3", expmap, "Exposure map (Level 3) from the CSC", move=True)
                data_package.store("fov3", fov, "FOV file (Level 3) from the CSC", move=True)

                # Make the data package read-only so we cannot change files by accident

                data_package.read_only = True

                # # Create directory named after obsid
                # if not os.path.exists(str(this_obsid)):
                #
                #     os.mkdir(str(this_obsid))
                #
                # # Move files in there
                #
                # os.rename(evtfile, os.path.join(str(this_obsid), evtfile))
                # os.rename(tsvfile, os.path.join(str(this_obsid), tsvfile))
                # os.rename(expmap, os.path.join(str(this_obsid), expmap))

            else:

                print "Region files do not exist for ObsID %s" % this_obsid
