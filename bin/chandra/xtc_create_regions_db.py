#!/usr/bin/env python

"""
Create the database of the regions, containing the R.A,Dec of each source and the path to its region file
"""

import argparse
import collections
import os
import re
import sys
import multiprocessing

import fitsio

from XtDac.ChandraUtils import find_files
from XtDac.ChandraUtils import work_within_directory
from XtDac.ChandraUtils.logging_system import get_logger
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename
from XtDac.ChandraUtils.configuration import get_configuration


def worker(region_file):

    header = fitsio.read_header(region_file, "SRCREG")

    ra = header['RA']
    dec = header['DEC']
    obsid = header['OBS_ID']

    assert ra is not None, "Cannot find R.A. in file %s" % region_file
    assert dec is not None, "Cannot find Dec. in file %s" % region_file
    assert obsid is not None, "Cannot find OBSID in file %s" % region_file

    try:

        ra = float(ra)
        dec = float(dec)

    except:

        raise RuntimeError("Cannot convert coordinates (%s,%s) to floats" % (ra, dec))

    try:

        obsid = int(obsid)

    except:

        raise RuntimeError("Cannot convert obsid %s to integer" % obsid)

    # Get the source name from the path

    try:

        g = re.match('(.+)/(CXOJ.+)/.+', region_file)

        name = g.groups()[1]

    except:

        raise RuntimeError("Cannot figure out the name of the source from the path %s" % (region_file))

    return (name, (ra, dec, obsid, os.path.relpath(region_file)))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create source database')

    parser.add_argument("-c", "--config_file", help="Path to the configuration file", type=str, required=True)

    # assumption = all level 3 region files and event file are already downloaded into same directory, the region_dir

    args = parser.parse_args()

    # Get logger
    logger = get_logger(os.path.basename(sys.argv[0]))

    # Get the configuration
    config = get_configuration(args.config_file)

    region_dir = sanitize_filename(config['region repository'])

    with work_within_directory.work_within_directory(region_dir):

        # Find all region files
        region_files = find_files.find_files('.', '*_reg3.fits.gz')

        logger.info("Found %s region files" % len(region_files))

        db = collections.OrderedDict()

        logger.info("Starting processing...")

        pool = multiprocessing.Pool(multiprocessing.cpu_count())

        try:

            for i, result in enumerate(pool.imap(worker, region_files)):

                sys.stderr.write("\r%i out of %i" % (i+1, len(region_files)))

                name, entry = result

                db[name] = entry

        except:

            raise

        finally:

            pool.close()
            sys.stderr.write("\n")

        logger.info("done")

        # Write the outfile

        with open('region_database.txt', 'w+') as f:

            f.write("#NAME RA DEC OBSID REGION_FILE\n")

            for source_name in db:
                ra, dec, obsid, region_file = db[source_name]

                f.write("%s %.5f %.5f %i %s\n" % (source_name, ra, dec, obsid, region_file))
