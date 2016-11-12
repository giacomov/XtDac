#!/usr/bin/env python

"""
Take each CCD file, already checked for hot pixels and cross matched for variable sources, and add remaining candidate
transients to master list
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils.data_package import DataPackage
from XtDac.ChandraUtils.run_command import CommandRunner
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Add filtered candidates to a masterlist of sources")

    parser.add_argument("--package", help="Data package containing the output of step 2",
                        required=True, type=str)
    parser.add_argument("--masterfile", help="Path to file containing list of transients in this set",
                        required=True, type=str)
    # parser.add_argument("--evtfile", help="Main event file for observation, used to get total exposure time",
    #                    required=True, type=str)


    # Get logger for this command

    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Instance the command runner

    runner = CommandRunner(logger)

    args = parser.parse_args()

    data_package = DataPackage(sanitize_filename(args.package))

    for bbfile_tag in data_package.find_all("ccd_?_check_var"):

        logger.info("Processing %s..." % bbfile_tag)

        bbfile = data_package.get(bbfile_tag).filename

        logger.info("(reading from file %s)" % bbfile)

        # get directory path and file name from input file arguments

        bb_file_path = sanitize_filename(bbfile)

        masterfile = sanitize_filename(args.masterfile)

        # evtfile = sanitize_filename(args.evtfile)

        # read BB data into array
        bb_data = np.array(np.recfromtxt(bb_file_path, names=True), ndmin=1)

        # number of rows of data
        bb_n = len(bb_data)

        if bb_n == 0:
            continue

        # column names
        existing_column_names = " ".join(bb_data.dtype.names)

        # with pyfits.open(evtfile) as evt:

        #    exposure = evt['EVENTS'].header['EXPOSURE']

        # If fresh list is to be started

        if not os.path.exists(masterfile):

            # new array of appropriate size
            master_data = np.empty([0, len(bb_data[0]) + 1])

        else:

            master_data = np.array(np.recfromtxt(masterfile, names=True), ndmin=1)

        with open('__temp.txt', "w") as temp:

            temp.write("# %s\n" % existing_column_names)

            i = 0

            while i < len(master_data):

                temp_list = []

                for j in xrange(1, len(bb_data.dtype.names) + 1):
                    temp_list.append(str(master_data[i][j]))

                line = " ".join(temp_list)

                temp.write("%s\n" % line)

                i += 1

            i = 0

            while i < len(bb_data):

                temp_list = []

                for j in xrange(len(bb_data.dtype.names)):
                    temp_list.append(str(bb_data[i][j]))

                line = " ".join(temp_list)

                temp.write("%s\n" % line)

                i += 1

        data_raw = np.array(np.recfromtxt("__temp.txt", names=True), ndmin=1)

        data_all = pd.DataFrame.from_records(data_raw)
        # Drop all duplicates

        data_unique = data_all.drop_duplicates(subset=['Candidate', 'Obsid', 'CCD'])

        # Keep only data satisfying the following conditions (empty at the moment)

        #idx = (data_unique['PSFfrac'] > 0)

        # Sort according to PSFfrac

        data_filtered = data_unique#[idx]

        master_data_sorted = sorted(data_filtered.values.tolist(), key=lambda data_row: data_row[-1], reverse=True)

        # Write data to text file

        with open(args.masterfile, "w") as f:

            f.write("# Rank %s\n" % existing_column_names)

            i = 0

            while i < len(master_data_sorted):

                temp_list = []

                temp_list.append(str(i + 1))

                for j in xrange(len(bb_data.dtype.names)):
                    temp_list.append(str(master_data_sorted[i][j]))

                line = " ".join(temp_list)

                f.write("%s\n" % line)

                i += 1

        os.remove('__temp.txt')
