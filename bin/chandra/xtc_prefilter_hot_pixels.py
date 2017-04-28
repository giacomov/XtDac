#!/usr/bin/env python

"""
Pre-filter event file for obvious hot pixels (to gain execution time). A further, deeper search for hot pixel will be
executed after the Bayesian Block stage
"""

import argparse
import os
import sys

import astropy.io.fits as pyfits
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import euclidean_distances

from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils.run_command import CommandRunner
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename

from XtDac.BayesianBlocks import bayesian_blocks


def neighbor_pixel_check(cluster_coords, max_distance=2):
    """
    Check that all events in the cluster have a maximum distance smaller than max_distance

    :param cluster_coords:
    :param max_distance:
    :return:
    """

    hot_pix_flag = True

    all_distances = euclidean_distances(cluster_coords, cluster_coords)
    if args.debug == "yes":
        print all_distances

    for distances in all_distances:

        for distance in distances:

            if distance < max_distance:
                pass

            else:
                hot_pix_flag = False
                break

        if hot_pix_flag == False:
            break

    return hot_pix_flag


def unique_rows(my_array):
    """
    Returns the data where all duplicated elements have been removed

    :param my_array: a n-dimensional array
    :return: unique entries in the n-dimensional array
    """

    # This is needed to avoid errors in the next
    my_array = np.ascontiguousarray(my_array)

    uniq = np.unique(my_array.view(my_array.dtype.descr * my_array.shape[1]))
    return uniq.view(my_array.dtype).reshape(-1, my_array.shape[1])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Check to see if transient candidates are actually hot pixels")

    parser.add_argument("--evtfile", help="Filtered CCD file", required=True, type=str)
    parser.add_argument("--outfile", help="Output file with events in hot pixels removed", required=True, type=str)
    parser.add_argument("--max_duration", help="Maximum duration to consider for bright pixels (default: 25)",
                        required=False, default=25, type=float)
    parser.add_argument("--debug", help="Debug mode? (yes or no)", required=False, default='no')

    # Get logger for this command

    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Instance the command runner

    runner = CommandRunner(logger)

    args = parser.parse_args()

    eventfile = sanitize_filename(args.evtfile)

    # Open event file

    tot_hot_pixels = 0

    with pyfits.open(eventfile, mode='update', memmap=False) as fits_file:

        # Get the data extension

        data = fits_file['EVENTS'].data

        # If this is a simulated data, there is no need for this processing
        if fits_file['EVENTS'].header.get("XTDACSIM"):

            logger.info("This is a simulated file. Skipping hot pixels preprocessing...")

        else:

            # Get unique CCD ids
            ccds = np.unique(fits_file['EVENTS'].data.ccd_id)

            # Get frame time

            tstart = fits_file['EVENTS'].header['TSTART']
            tstop = fits_file['EVENTS'].header['TSTOP']

            for ccd in ccds:

                logger.info("Processing CCD %s" % ccd)

                # Select all events in this CCD
                ccd_idx = (data.ccd_id == ccd)

                chipx = data.chipx[ccd_idx]
                chipy = data.chipy[ccd_idx]
                time = data.time[ccd_idx]

                coords = np.vstack([chipx, chipy]).T

                ucoords = unique_rows(coords)

                # For each non-empty pixel do a bayesian block analysis
                for i, (cx, cy) in enumerate(ucoords):

                    if (i+1) % 10000 == 0:

                        logger.info("%s out of %s" % (i+1, coords.shape[0]))

                    time_stamps = data.time[ccd_idx & (data.chipx == cx) & (data.chipy == cy)]

                    # Do not try if there are only 5 events in the whole observation in this pixel

                    if time_stamps.shape[0] < 5:

                        continue

                    blocks = bayesian_blocks(time_stamps, tstart, tstop, 1e-3)

                    if len(blocks) > 2:

                        for t1,t2 in zip(blocks[:-1], blocks[1:]):

                            duration = t2 - t1

                            if duration < args.max_duration:

                                # Check if this is a bright pixel

                                # Select all events within this time interval
                                time_idx = (time >= t1) & (time <= t2)

                                # Select all events in this time interval in the surrounding pixels

                                d = np.sqrt((chipx-float(cx))**2 + (chipy-float(cy))**2)

                                this_idx = d == 0

                                neighbor_idx = time_idx & (d < 2) & ~this_idx

                                if np.sum(neighbor_idx) == 0:

                                    # No events in surrounding pixels. This is likely a hot pixel
                                    logger.info(" @ (%s, %s), interval: %1.f - %.1f s "
                                                "(%.1f s, %i evts)" % (cx, cy, t1, t2, duration, np.sum(this_idx)))

                                    # Flag the events
                                    events_idx = ccd_idx & (data.chipx == cx) & (data.chipy == cy) & \
                                                 (data.time >= t1) & (data.time <= t2)

                                    data.pha[events_idx] = -1

                                    tot_hot_pixels += 1


    # # Count how many events we are filtering out
    n_filtered = np.sum(data.pha == -1)
    #
    n_tot = int(data.shape[0])

    logger.info("Found %s non-unique hot pixels" % tot_hot_pixels)
    logger.info(
        "Filtering out %s events out of %s (%.2f percent)" % (n_filtered, n_tot, float(n_filtered) / n_tot * 100.0))

    outfile = sanitize_filename(args.outfile)

    cmd_line = "ftcopy %s[EVENTS][pha!=-1] %s  copyall=true clobber=yes" % (eventfile, outfile)

    runner.run(cmd_line)
