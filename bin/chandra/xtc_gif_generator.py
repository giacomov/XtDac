#!/usr/bin/env python

"""
Generate lightcurves for each candidate given a list of candidates
"""

# Set to a non-interactive matplotlib backend

import matplotlib
matplotlib.use("agg")


import argparse
import os
import sys
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import seaborn as sbs
from matplotlib.colors import LogNorm
from astropy.convolution import convolve_fft, convolve, Gaussian2DKernel
from matplotlib.animation import ArtistAnimation, FFMpegWriter

from XtDac.ChandraUtils.find_files import find_files
from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils.run_command import CommandRunner
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename

from XtDac.DivideAndConquer import XMMWCS
from XtDac.DivideAndConquer.HardwareUnit import hardwareUnitFactory

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Generate gifs to visualize transient')

    parser.add_argument("--masterfile", help="Path to file containing list of transients",
                        required=True, type=str)
    parser.add_argument("--data_path", help="Path to directory containing data", required=False, type=str, default='.')

    parser.add_argument("--verbose-debug", action='store_true')
    parser.add_argument("--cleanup", action='store_true')

    # Get the logger
    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Get the command runner
    runner = CommandRunner(logger)

    args = parser.parse_args()

    data_path = sanitize_filename(args.data_path)
    masterfile = sanitize_filename(args.masterfile)

    transient_data = np.array(np.recfromtxt(masterfile, names=True), ndmin=1)

    for transient in transient_data:

        obsid = transient['Obsid']
        ccd = transient['CCD']
        candidate = transient['Candidate']
        tstart = transient['Tstart']
        tstop = transient['Tstop']

        duration = tstop-tstart

        event_files = find_files(data_path, "ccd_%s_*_filtered_*.fits" % (ccd))

        assert len(event_files)==1, "Couldn't find event file. " \
                                    "I was looking for %s" % ("ccd_%s_*_filtered_*.fits" % (ccd))

        event_file = event_files[0]

        # get start and stop time of observation
        with pyfits.open(event_file, memmap=False) as fits_file:

            # Get the start of the first GTI and the stop of the last one
            gti_starts = []
            gti_stops = []
            for ext in fits_file[1:]:

                if ext.header['EXTNAME'] == 'GTI':

                    gti_starts.append(ext.data.field("START").min())
                    gti_stops.append(ext.data.field("STOP").max())

            frame_time = fits_file['EVENTS'].header['TIMEDEL']

            tmin = min(gti_starts)
            tmax = max(gti_stops)

            # Get minimum and maximum X and Y, so we use always the same binning for the images
            xmin, xmax = fits_file['EVENTS'].data.field("X").min(), fits_file['EVENTS'].data.field("X").max()
            ymin, ymax = fits_file['EVENTS'].data.field("Y").min(), fits_file['EVENTS'].data.field("Y").max()

            # Get position and transform to X,Y

            ra = transient['RA']
            dec = transient['Dec']

            wcs = XMMWCS.XMMWCS(event_file, fits_file['EVENTS'].data.field("X"), fits_file['EVENTS'].data.field("Y"))
            transient_x, transient_y = wcs.sky2xy([[ra, dec]])[0]

        hwu = hardwareUnitFactory(event_file)

        logger.info("Duration: %s" %duration)
        logger.info("Tmin: %s" % tmin)
        logger.info("Tmax: %s" %tmax)
        logger.info("frame time: %s" % frame_time)
        logger.info("Obs Time: %s" %(tmax-tmin))
        logger.info("X,Y = (%.3f, %.3f)" % (transient_x, transient_y))

        # Use the interval before the transient and the interval after the transient, as well as of course
        # the interval of the transient itself

        intervals = [tmin]

        # Add tstart only if it is sufficiently different from tmin (so we don't get an empty interval when the
        # transient is right at the beginning)

        if abs(tstart - tmin) > (frame_time + frame_time / 2.0):

            intervals.append(tstart)

        intervals.append(tstop)

        # Add tmax only if it is sufficiently different from tstop, so we don't get an empty interval when the
        # transient is right at the end)

        if abs(tmax - tstop) > (frame_time + frame_time / 2.0):

            intervals.append(tmax)

        intervals = sorted(intervals)

        deltas = np.array(intervals[1:]) - np.array(intervals[:-1])

        evt_name, evt_file_ext = os.path.splitext(os.path.basename(event_file))

        # Create individual time interval files
        images = []

        for i in range(len(intervals)-1):

            outfile = "cand_%s_%s_TI_%s%s" %(candidate, evt_name, i+1, evt_file_ext)

            cmd_line = 'ftcopy \"%s[(TIME >= %s) && (TIME <= %s)]\" %s clobber=yes ' \
                       % (event_file, intervals[i], intervals[i+1], outfile)

            runner.run(cmd_line)

            images.append(outfile)

        #create a list of frames that will be animated into a gif
        frames = []

        # Prepare bins
        # xbins = np.linspace(xmin, xmax, 300)
        # ybins = np.linspace(ymin, ymax, 300)

        conv_kernel_size = max(np.ceil(float(transient['PSF_sizearcsec']) / hwu.getPixelScale()), 5)

        # Create bins of size 1 (i.e., one pixel per bin)

        xbins = np.arange(transient_x - 5 * conv_kernel_size, transient_x + 5 * conv_kernel_size + 1, 1)
        ybins = np.arange(transient_y - 5 * conv_kernel_size, transient_y + 5 * conv_kernel_size + 1, 1)

        logger.info("Image will be %i x %i pixels" % (xbins.shape[0], ybins.shape[0]))

        fig = plt.figure()

        sub = fig.add_subplot(111)
        sub.set_title("ObsID %s, CCD %s \nTstart = %s, Tstop = %s (%i events)" % (obsid, ccd, tstart, tstop,
                                                                                  float(transient['N_events'])))

        for i, image in enumerate(images):

            # Compute interval duration
            dt = intervals[i + 1] - intervals[i]

            #get x and y coordinates of image from fits file
            data = pyfits.getdata(image)
            x = data.field("x")
            y = data.field("y")

            #create 2D histogram data from it
            sbs.set(rc={'axes.facecolor': 'black', 'axes.grid': False})

            hh, X, Y = np.histogram2d(x, y, bins=[xbins, ybins])

            #smooth data
            gauss_kernel = Gaussian2DKernel(stddev=max(conv_kernel_size / 8.0, 2.0))
            smoothed_data_gauss = convolve_fft(hh, gauss_kernel, normalize_kernel=True)

            if x.shape[0] > 0:

                img = sub.imshow(smoothed_data_gauss, cmap='hot', animated=True, origin='lower')

            else:

                # No events in this image. Generate an empty image
                img = sub.imshow(smoothed_data_gauss, cmap='hot', animated=True, origin='lower')

            # Draw PSF circle
            circ = Circle((smoothed_data_gauss.shape[0]/2.0 + 1, smoothed_data_gauss.shape[1]/2.0 + 1),
                          float(transient['PSF_sizearcsec']) / hwu.getPixelScale(),
                          fill=False, lw=2, color='green')
            sub.add_patch(circ)

            text = sub.annotate("%i" % (i+1), xy=(0.5, 0.03),
                                xycoords='figure fraction',
                                annotation_clip=False)

            text3 = sub.annotate("Duration: %.2f s" % (dt), xy=(0.55, 0.13),
                                 xycoords='figure fraction',
                                 annotation_clip=False, color='green')

            sub.set_yticklabels([])
            sub.set_xticklabels([])

            #add this image to list of frames
            frames.append([img, text, text3, circ])

            # Remove the image
            if args.cleanup:

                os.remove(image)

        #animate and save gif
        logger.info("Creating gif ObsID %s, CCD %s, Candidate %s...\n" %(obsid, ccd, candidate))
        anim = ArtistAnimation(fig, frames, interval=2000)
        anim.save("%s_cand_%s.gif" %(evt_name, candidate), writer='imagemagick')

        plt.close()


