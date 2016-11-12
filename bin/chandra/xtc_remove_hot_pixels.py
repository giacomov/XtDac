#!/usr/bin/env python

"""
Check to see if candidate transients are actually hot pixels
"""

import argparse
import glob
import os
import sys

import astropy.io.fits as pyfits
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import euclidean_distances

from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils.run_command import CommandRunner

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Check to see if transient candidates are actually hot pixels")

    parser.add_argument("--obsid", help="Observation ID Numbers", type=int, required=True)
    parser.add_argument("--evtfile", help="Filtered CCD file", required=True)
    parser.add_argument("--bbfile", help="Text file for one CCD with transient candidates listed", required=True)
    parser.add_argument("--outfile", help="Name of output text file file", required=True)
    parser.add_argument("--debug", help="Debug mode? (yes or no)", required=True)

    # Get logger for this command

    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Instance the command runner

    runner = CommandRunner(logger)

    args = parser.parse_args()

    # Find region files for each candidate transient

    evtfile = os.path.abspath(os.path.expandvars(os.path.expanduser(args.evtfile)))
    bbfile = os.path.abspath(os.path.expandvars(os.path.expanduser(args.bbfile)))

    # read BB data into array
    bb_data = np.array(np.recfromtxt(bbfile, names=True), ndmin=1)

    evt_file_name = os.path.splitext(os.path.basename(evtfile))[0]

    reg_files = glob.glob('%s_candidate*reg' % evt_file_name)


    # make sure files are sorted

    def extract_number(s):
        return int(os.path.splitext(s)[0].split("_")[-1])


    reg_files_sorted = sorted(reg_files, key=extract_number)

    # Find CCD number based on file name

    names = os.path.splitext(evt_file_name)[0].split("_")
    idx = names.index("ccd")

    ccd_num = names[idx + 1]


    def neighbor_pixel_check(cluster_coords):

        hot_pix_flag = True

        all_distances = euclidean_distances(cluster_coords, cluster_coords)
        if args.debug == "yes":
            print all_distances

        for distances in all_distances:

            for distance in distances:

                if distance < 2:
                    pass

                else:
                    hot_pix_flag = False
                    break

            if hot_pix_flag == False:
                break

        return hot_pix_flag


    with open(args.outfile, "w") as f:

        # Pre-existing column names
        existing_column_names = " ".join(bb_data.dtype.names)

        f.write("# Candidate Obsid CCD %s Duration N_events Hot_Pixel_Flag\n" % existing_column_names)

        for n, reg_file in enumerate(reg_files_sorted):

            tstart = bb_data['Tstart'][n]
            tstop = bb_data['Tstop'][n]

            temp_reg_file = "temp_reg_%s.fits" % (n + 1)

            # Make temporary fits region file with region determined by xtdac
            cmd_line = 'ftcopy \"%s[EVENTS][regfilter(\'%s\') && (TIME >= %s) && (TIME <= %s)]\" %s clobber=yes ' \
                       % (evtfile, reg_file, tstart, tstop, temp_reg_file)

            runner.run(cmd_line)

            # Initialize hot_pix_flag to False
            hot_pix_flag = False

            # Open temporary region file to get chip coordinates
            with pyfits.open(temp_reg_file, memmap=False) as reg:
                chipx = reg['EVENTS'].data.chipx
                chipy = reg['EVENTS'].data.chipy

            # Count number of events in the region in this time interval
            n_events = chipx.shape[0]

            if len(chipx) == 0:

                print "\n\n%s had no events between %s and %s!\n\n" % (reg_file, tstart, tstop)

            elif len(chipx) < 15:

                # Organize into list of coordinate pairs
                coords = []

                for i in range(len(chipx)):
                    temp = []
                    temp.append(chipx[i])
                    temp.append(chipy[i])
                    coords.append(temp)

                # Run DBSCAN
                db = DBSCAN(eps=2, min_samples=2).fit(coords)
                l = db.labels_
                csi = db.core_sample_indices_
                c = db.components_

                if args.debug == "yes":
                    print coords
                    print l

                # Run this if more than one cluster found
                if len(np.unique(l)) > 1:

                    labels_unique = np.unique(l)

                    if args.debug == "yes":
                        print labels_unique

                    if labels_unique[0] == -1:
                        labels_unique = labels_unique[1::]

                    if args.debug == "yes":
                        print labels_unique

                    hot_pix_flags = []
                    org_data_idx = []

                    # Create a list of indices corresponding to each label
                    # [[indices of items in cluster 0],[indices of items in cluster 1],...]

                    for i in labels_unique:

                        temp = []
                        for a, b in enumerate(l):

                            if b == i:
                                temp.append(a)

                        org_data_idx.append(temp)

                    # Get x and y coordinates for each label
                    for i in range(len(org_data_idx)):
                        x = []
                        y = []
                        cluster_coords = []

                        for idx in org_data_idx[i]:
                            cluster_coords.append(coords[idx])
                            x.append(coords[idx][0])
                            y.append(coords[idx][1])

                        if args.debug == "yes":
                            print cluster_coords

                        npx = np.array(x)
                        npy = np.array(y)

                        npxl = len(np.unique(npx))
                        npyl = len(np.unique(npy))

                        # Check if coordinates are all the same and adjust hot_pix_flag accordingly

                        if npxl == 1 and npyl == 1:
                            hot_pix_flags.append(True)

                        # If not, check for neighboring coordinates
                        else:

                            hot_pix_flag = neighbor_pixel_check(cluster_coords)

                            hot_pix_flags.append(hot_pix_flag)

                    # If every label is flagged as a hot pixel, flag transient as hot pixel
                    hot_pix_flags = np.array(hot_pix_flags)

                    flags = np.unique(hot_pix_flags)

                    if len(flags) == 1 and flags[0] == True:
                        hot_pix_flag = True

                    if args.debug == "yes":
                        print hot_pix_flags
                        print hot_pix_flag

                # If only one cluster is found

                else:
                    x = []
                    y = []
                    for i in range(len(coords)):
                        x.append(coords[i][0])
                        y.append(coords[i][1])

                    npxl = len(np.unique(np.array(x)))
                    npyl = len(np.unique(np.array(y)))

                    # Check that all coordinates in cluster are the same
                    if npxl == 1 and npyl == 1:
                        hot_pix_flag = True

                    # If all coordinates are not same, check for neighboring pixels

                    if hot_pix_flag == False:
                        hot_pix_flag = neighbor_pixel_check(coords)

                    if args.debug == "yes":
                        print hot_pix_flag

            # Write to outfile

            temp_list = []

            temp_list.append(str(n + 1))
            temp_list.append(str(args.obsid))
            temp_list.append(str(ccd_num))

            for j in range(len(bb_data.dtype.names)):
                temp_list.append(str(bb_data[n][j]))

            # Fill Duration column
            temp_list.append("%.1f" % (tstop - tstart))

            # Fill N_events column
            temp_list.append("%i" % (n_events))

            # Fill "Hot_Pixel_Flag" column

            temp_list.append(str(hot_pix_flag))

            line = " ".join(temp_list)

            f.write("%s\n" % line)

            if args.debug == "no":
                os.remove(temp_reg_file)

    if args.debug == "yes":
        print "NOTE: Debug mode, temporary region files not deleted"
