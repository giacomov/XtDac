#!/usr/bin/env python

"""
Take each CCD file, already checked for hot pixels and cross matched for variable sources, and add remaining candidate
transients to master list
"""

import argparse
import os
import sys

import astropy.io.fits as pyfits
from astropy import wcs as pywcs
import astropy.table

import numpy as np

from XtDac.ChandraUtils.run_command import CommandRunner
from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Simulate a background-only observation")

    parser.add_argument("--bkgmap", help="Background map file",
                        required=True, type=str)
    parser.add_argument("--evtfile", help="Filtered event file with sources removed",
                        required=True, type=str)
    parser.add_argument("--asolfile", help="Aspect solution file",
                        required=True, type=str)
    parser.add_argument("--outfile", help="Output file containing the simulated events",
                        required=True, type=str)

    # Get logger for this command

    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Parse arguments

    args = parser.parse_args()

    # Sanitize file names and make sure they exist

    bkgfile = sanitize_filename(args.bkgmap)
    asolfile = sanitize_filename(args.asolfile)
    evtfile = sanitize_filename(args.evtfile)
    outfile = sanitize_filename(args.outfile)

    for filename in [bkgfile, asolfile, evtfile]:

        assert os.path.exists(filename), "File %s does not exist" % filename

    logger.info("Reading background image %s..." % bkgfile)

    data, header = pyfits.getdata(bkgfile, 0, header=True)

    # Build WCS to convert back and forth from image coordinates to X,Y coordinates

    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = [header['CRPIX1P'], header['CRPIX2P']]
    wcs.wcs.cdelt = [header['CDELT1P'], header['CDELT2P']]
    wcs.wcs.crval = [header['CRVAL1P'], header['CRVAL2P']]
    wcs.wcs.ctype = [header['CTYPE1P'], header['CTYPE2P']]

    # Generate a new background image with Poisson noise

    logger.info("Generating simulated dataset...")

    new_data = np.random.poisson(data)

    # Get start and stop of the observation, and its duration
    tstart = header['TSTART']
    tstop = header['TSTOP']
    duration = tstop - tstart

    # Now unravel the data creating an unbinned event list

    xs = []
    ys = []
    ts = []

    for i in range(data.shape[0]):

        for j in range(data.shape[1]):

            # NOTE: most of the pixels in the new_data image area actually 0
            # (there is a large padding around the image)

            n_events_to_generate = new_data[i, j]

            if n_events_to_generate > 0:

                # Convert from image coordinates to physical coordinates (i.e., X and Y)

                pixcrd = np.array([[i, j]], np.float_)

                xy = wcs.all_pix2world(pixcrd, 1)[0]

                # Generate time of the events according to a uniform Poisson distribution
                # with rate given by the background map for this pixel

                rate = data[i, j] / duration

                new_times = np.cumsum(np.random.exponential(1.0 / rate, n_events_to_generate)) + tstart

                # Append the new times and the new X,Y coordinates

                ts.extend(new_times)

                # Note the swap in xy (1 is x and 0 is y). This is not an error. It is required to
                # match the WCS in the event file

                xs.extend([xy[1]] * n_events_to_generate)
                ys.extend([xy[0]] * n_events_to_generate)

    # Transform into arrays

    ts = np.array(ts)
    xs = np.array(xs)
    ys = np.array(ys)

    assert np.sum(new_data) == xs.shape[0]

    # Make sure we stay within tstart and tstop

    idx = (ts >= tstart) & (ts <= tstop)
    ts = ts[idx]
    xs = xs[idx]
    ys = ys[idx]

    # Now read in the event file and substitute the x,y and time values of a subsample of events
    # This way we will keep a realistic spectrum for the background

    logger.info("Generating FITS file...")

    new_table = astropy.table.Table.read(evtfile, 'EVENTS')

    # Get the list of CCDs present
    unique_ccds = np.unique(new_table['ccd_id'])

    # Pick up randomly events in the file, we will overwrite their x, y and time
    # NOTE: randint will happily reuse events (which is not a problem) so we can obtain more events than
    # the one contained in the event file

    n_sim_evt = len(xs)

    idx = np.random.randint(0, len(new_table), n_sim_evt)

    new_table = new_table[idx]

    # Now overwrite their coordinates and times with the new ones

    new_table['x'] = xs
    new_table['y'] = ys
    new_table['time'] = ts

    # Now write them in a temporary file

    # Keep track of temp files
    temporary_files = []

    temp_evt_file = 'test.fits'

    new_table.write(temp_evt_file, overwrite=True)

    temporary_files.append(temp_evt_file)

    # Now divide in CCDs

    runner = CommandRunner(logger)

    # Generate region file containing the chips

    region_file = '__chips_reg.fits'
    cmd_line = 'skyfov %s %s aspect=%s clobber=yes' % (evtfile, region_file, asolfile)

    runner.run(cmd_line)

    temporary_files.append(region_file)

    # Now divide all chips. Note that this will count some events twice (at the border of the CCDs),
    # as there is overlap between the CCD regions in sky coordinates due to small motions in the pointing
    # This is not a problem, as we will "uniquefy" them at the end

    ccd_file_list = "__ccd_files.txt"

    with open(ccd_file_list, "w+") as f:

        for ccd_id in unique_ccds:

            ccd_file = '__ccd_%i.fits' % ccd_id

            cmd_line = 'dmcopy "%s[sky=region(%s[ccd_id=%i])]" %s clobber=yes' % (temp_evt_file, region_file,
                                                                                  ccd_id, ccd_file)

            runner.run(cmd_line)

            # Now open the file and update the ccd_id column
            with pyfits.open(ccd_file, mode='update') as fitsf:
                fitsf['EVENTS'].data.ccd_id[:] = ccd_id

            f.write("%s\n" % ccd_file)

            temporary_files.append(ccd_file)

    temporary_files.append(ccd_file_list)

    # Merge them back
    cmd_line = "ftmerge @%s %s clobber=yes copyall=yes" % (ccd_file_list, temp_evt_file)

    runner.run(cmd_line)

    # Sort and remove duplicates

    cmd_line = "ftsort %s[EVENTS] %s 'time' unique=yes clobber=yes copyall=yes memory=yes" % (temp_evt_file,
                                                                                              outfile)

    runner.run(cmd_line)

    # Remove temporary files
    for filename in temporary_files:

        os.remove(filename)