#!/usr/bin/env python

"""
Download files by observation ID
"""

import argparse
import os
import sys
import multiprocessing

from XtDac.ChandraUtils import find_files
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename
from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils.data_package import DataPackage
from XtDac.ChandraUtils.run_command import CommandRunner
from XtDac.ChandraUtils.work_within_directory import work_within_directory
from XtDac.ChandraUtils.configuration import get_configuration


def worker(this_obsid):

    regdir_this_obsid = os.path.join(region_repository, str(this_obsid))

    if os.path.exists(regdir_this_obsid):

        # This could download more than one observation segment for this obsid

        cmd_line = "xtc_download_by_obsid.py --obsid %d" % this_obsid

        runner.run(cmd_line)

        # Get the downloaded files
        evtfiles = find_files.find_files(os.getcwd(), 'acisf%05i*evt3.fits' % int(this_obsid))

        logger.info("Found %s event files" % len(evtfiles))

        for evtfile in evtfiles:
            # Get the root of the evt3 file name
            # The evt3 file name is like acisf01578_001N001_evt3.fits.gz,
            # where the 001 is the observation segment
            name_root = "_".join(os.path.basename(evtfile).split("_")[:-1])  # this is "acisf01578_001N001"
            obsid_identifier = name_root.replace("acisf", "").split("N")[0]  # this is 01578_001

            logger.info("Processing %s" % obsid_identifier)

            # Find exposure map and fov file
            expmaps = find_files.find_files(os.getcwd(), "%s*exp3.fits*" % name_root)
            fovs = find_files.find_files(os.getcwd(), "%s*fov3.fits*" % name_root)
            tsvfiles = find_files.find_files(os.getcwd(), "%s.tsv" % obsid_identifier)
            bkgmaps = find_files.find_files(os.getcwd(), "*%s*bkgimg3.fits*" % this_obsid)
            asol_files = find_files.find_files(os.getcwd(), '*asol*.fits.gz')

            assert len(expmaps) == 1, "Wrong number of exposure maps for event file %s" % evtfile
            assert len(fovs) == 1, "Wrong number of fov files for event file %s" % evtfile
            assert len(tsvfiles) == 1, "Wrong number of tsv files for obsid %s" % this_obsid
            assert len(bkgmaps) == 1, "Wrong number of bkg files for obsid %s" % this_obsid
            assert len(asol_files) == 1, "Wrong number of asol files for obsid %s" % this_obsid


            tsvfile = tsvfiles[0]
            expmap = expmaps[0]
            fov = fovs[0]
            bkgmap = bkgmaps[0]
            asol = asol_files[0]

            logger.info("Found tsv file: %s" % tsvfile)
            logger.info("Found expmap: %s" % expmap)
            logger.info("Found fov file: %s" % fov)
            logger.info("Found bkg map file: %s" % bkgmap)
            logger.info("Found asol file: %s" % asol)

            logger.info("Creating data package %s" % obsid_identifier)

            data_package = DataPackage(obsid_identifier, create=True)

            data_package.store("evt3", evtfile, "Event file (Level 3) from the CSC", move=True)
            data_package.store("tsv", tsvfile, "TSV file from the CSC", move=True)
            data_package.store("exp3", expmap, "Exposure map (Level 3) from the CSC", move=True)
            data_package.store("fov3", fov, "FOV file (Level 3) from the CSC", move=True)
            data_package.store("bkgmap", bkgmap, "Background map (Level 3) from the CSC", move=True)
            data_package.store("asol", asol, "Aspect solution file from the CSC", move=True)

            logger.info("done")

            # Make the data package read-only so we cannot change files by accident
            logger.info("Making it read-only")
            data_package.read_only = True

    else:

        logger.error("Region files do not exist for ObsID %s" % this_obsid)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run following steps: download files by obsid')

    parser.add_argument("-c", "--config_file", help="Path to the configuration file", type=str, required=True)

    parser.add_argument("-o", "--obsid", help="Observation ID Numbers", type=int, required=True, nargs='+')

    parser.add_argument("-n", "--n_processes", help="Number of parallel processes to use",
                        type=int, required=False, default=1)

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

        pool = multiprocessing.Pool(args.n_processes)

        try:

            for i, _ in enumerate(pool.imap(worker, args.obsid)):

                sys.stderr.write("\n\n%i out of %s\n\n" % (i+1, len(args.obsid)))

        except:

            raise

        finally:

            pool.close()
