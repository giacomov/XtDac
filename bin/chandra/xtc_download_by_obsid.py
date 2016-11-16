#!/usr/bin/env python

"""
Download Level 3 Event File (evt3), exposure maps and FOV file given inputs of Observation IDs.

It also applies the r4_header_update script to the event file so that it is updated to the
Reprocessing 4 version of data.

Make sure CIAO is running before running this script
"""

import argparse
import os
import shutil
import sys
import warnings

from XtDac.ChandraUtils.work_within_directory import work_within_directory
from XtDac.ChandraUtils import find_files
from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils.run_command import CommandRunner
from XtDac.ChandraUtils.configuration import get_configuration

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Download event files and exposure map from the Chandra catalog and '
                                                 'place them in the data repository defined in the configuration file')

    parser.add_argument("-o", "--obsid", help="Observation ID Numbers (one or more, separated by a space)",
                        type=int, required=True, nargs="+")

    # Some of the commands need a temporary directory to write files, defined in ASCDS_WORK_PATH.
    # Enforce that such a variable is defined in the current environment
    if os.environ.get("ASCDS_WORK_PATH") is None:
        raise RuntimeError("You need to set the env. variable ASCDS_WORK_PATH to a writeable temporary directory")

    # Get logger for this command

    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Instance the command runner

    runner = CommandRunner(logger)

    # default directory = current one

    args = parser.parse_args()

    for obsid in args.obsid:

        # Create temporary directory

        temp_dir = "__temp"

        if not os.path.exists(temp_dir):

            os.makedirs("__temp")

        else:

            warnings.warn("The directory %s already exists" % temp_dir)

        # Move there and download files

        with work_within_directory(temp_dir):

            work_dir = os.path.abspath(os.getcwd())

            # Download exposure map

            cmd_line = ("obsid_search_csc obsid=%d download=all outfile=%d.tsv filetype=exp,evt,fov "
                        "mode=h clobber=yes verbose=0 "
                        "columns=m.ra,m.dec,o.theta,m.extent_flag,m.var_flag"
                        % (obsid, obsid))

            runner.run(cmd_line)

            # Download ancillary files needed by the r4_header_update script
            cmd_line = "download_chandra_obsid %d asol,pbk -q" % (obsid)

            runner.run(cmd_line)

        # get paths of files
        evt3_files = find_files.find_files(work_dir, '*%s*evt3.fits.gz' % obsid)
        tsv_files = find_files.find_files(work_dir, '*%s*.tsv' % obsid)

        if len(evt3_files) == 0:

            raise RuntimeError("Could not find any event file. Maybe download failed?")

        else:

            # The asol and pbk files contain the start time in their name, as in
            # acisf102646229N003_pbk0.fits.gz, thus we find all of them here and we order
            # them

            asol_files = find_files.find_files(work_dir, '*asol*.fits.gz')
            pbk_files = find_files.find_files(work_dir, '*pbk*.fits.gz')

            # Order them
            def orderer(name):

                return int(os.path.basename(name).split("N")[0][5:])

            asol_files = sorted(asol_files, key=orderer)
            pbk_files = sorted(pbk_files, key=orderer)

            assert len(asol_files) > 0 and len(pbk_files) > 0, "Couldn't find ASOL or PBK files for obsid %s" % obsid

            logger.info("Found %s observation segments" % len(evt3_files))

            for evt3, pbk, asol in zip(evt3_files, pbk_files, asol_files):

                # Get the root of the evt3 file name
                # The evt3 file name is like acisf01578_001N001_evt3.fits.gz,
                # where the 001 is the observation segment
                name_root = "_".join(os.path.basename(evt3).split("_")[:-1])  # this is "acisf01578_001N001"

                # Find exposure map and fov file
                exp_files = find_files.find_files(work_dir, "%s_*exp3.fits.gz" % name_root)
                fov_files = find_files.find_files(work_dir, "%s_*fov3.fits.gz" % name_root)

                assert (len(exp_files) > 0) and len(fov_files) > 0, "Couldn't find exposure map or fov file " \
                                                                    "corresponding to event file %s" % evt3

                assert len(exp_files) == 1 and len(fov_files) == 1, "More than one exposure map or fov file for event " \
                                                                    "file %s" % evt3

                exp = exp_files[0]
                fov = fov_files[0]

                # The r4_header_update script cannot run on a compressed fits file, so decompress the eventfile

                cmd_line = "gunzip %s" % evt3

                runner.run(cmd_line)

                evt3 = evt3.replace("fits.gz", "fits")

                # Run reprocessing
                cmd_line = "r4_header_update infile=%s pbkfile=%s asolfile=%s" % (evt3, pbk, asol)

                runner.run(cmd_line)

                # Copy tsv file
                obs_identifier = name_root.replace("acisf","").split("N")[0]  # this is 01578_001

                new_tsv = "%s.tsv" % obs_identifier

                shutil.copy2(tsv_files[0], new_tsv)

                # move files

                for this_file in [evt3, exp, fov]:

                    os.rename(this_file, os.path.basename(this_file))

            # remove the tsv file
            os.remove(tsv_files[0])

        # Remove temp directory

        shutil.rmtree(temp_dir)
