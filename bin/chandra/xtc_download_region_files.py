#!/usr/bin/env python

import subprocess
import numpy as np
import os
import sys
import argparse
import multiprocessing
import time
import glob

from XtDac.data_files import get_data_file_path
from XtDac.ChandraUtils.work_within_directory import work_within_directory
from XtDac.ChandraUtils.logging_system import get_logger
from XtDac.ChandraUtils.configuration import get_configuration
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename
from XtDac.ChandraUtils.find_files import find_files


def worker(obsid):

    obsid = int(obsid)

    cmd_line = 'obsid_search_csc obsid=%s outfile=%s.tsv download=all ' \
               'filetype=reg mode=h clobber=yes verbose=0' % (obsid, obsid)

    for j in range(10):

        try:

            subprocess.check_call(cmd_line, shell=True)

        except subprocess.CalledProcessError:

            # Wait 60 seconds then retry

            time.sleep(60.0)

            continue

        else:

            # Check that we have all files
            with open("%s.tsv" % obsid) as f:

                lines = f.readlines()

            # Remove comments
            lines = filter(lambda x:x[0]!='#', lines)

            # First line after the comment is the column head, so we need to subtract one
            number_of_sources = len(lines) - 1

            number_of_files = len(find_files("%s" % obsid, "*_reg3.fits.gz"))

            if number_of_sources != number_of_files:

                logger.error("Number of downloaded files different from number of sources. Retrying...")

                continue

            else:

                # All good, go to next obsid

                break

    if not os.path.exists(str(obsid)):


        return None

    else:

        return obsid

if __name__=="__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description='Download region files from the Chandra Source Catalog in the region '
                                                 'repository')

    parser.add_argument("-c", "--config_file", help="Path to the configuration file", type=str, required=True)

    args = parser.parse_args()

    # Get logger
    logger = get_logger(os.path.basename(sys.argv[0]))

    # Get the configuration
    config = get_configuration(args.config_file)

    # Get region repository path
    region_repository = sanitize_filename(config['region repository'])

    logger.info("Saving regions to %s" % region_repository)

    # Read observation ids

    obs_list_file = get_data_file_path("obs_in_csc1.1.txt")

    logger.info("Reading list of observations from %s" % obs_list_file)

    obsids = np.genfromtxt(obs_list_file)

    with work_within_directory(region_repository, create=True):

        pool = multiprocessing.Pool(multiprocessing.cpu_count())

        good = []
        bad = []

        try:

            for i, res in enumerate(pool.imap(worker, obsids)):

                if res is not None:

                    good.append(res)

                else:

                    bad.append(obsids['obsid'][i])

                logger.info("%s out of %s (%s bad)" % (i + 1, len(obsids), len(bad)))

        except:

            raise

        finally:

            pool.close()

    logger.info("Downloaded region files for %s observations out of %s" % (len(good), len(obsids)))

    if len(bad) > 0:

        logger.error("\n\nWARNING: could not download region files for %s observations!" % len(bad))
        logger.error("Specifically:")
        logger.error(",".join(map(lambda x:str(int(x)), bad)))


