#!/usr/bin/env python

# Author: Giacomo Vianello (giacomov@stanford.edu)

##############
# Fixed Bin Search for transients
##############

# Select backend before anything else,
# otherwise the other import will override this
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

# Set up the logger (its verbosity level will be changed later on
# using the value provided by the user)
import logging

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)-15s %(message)s")

log = logging.getLogger("main")

import argparse
import os, sys
import numpy
import warnings

try:

    import astropy.io.fits as pyfits

except:

    import pyfits

from XtDac.bin.xtdac import fitsRegions
from XtDac.bin.xtdac import FixedBinSearch
from XtDac.bin.xtdac import XMMWCS
from XtDac.bin.xtdac import HardwareUnit

import scipy.stats


def validFITSfile(arg):
    if not os.path.isfile(arg):
        log.error("The file %s does not exist!" % arg)
        sys.exit(-1)

    # Try to open it to see if it is a FITS file
    try:
        pyfits.open(arg)
    except:
        log.error("The file %s is not a valid FITS file" % arg)
    return arg


parser = argparse.ArgumentParser(description="Fixed Bin Search for transients")

parser.add_argument("-t", "--dt",
                    help="Size of the fixed time bin",
                    required=True, type=float)

parser.add_argument("-e", "--eventfile",
                    help="Event file to use (already cleaned and selected in" +
                         " energy and quadrant/CCD), and containing the region definitions",
                    required=True, type=validFITSfile)

parser.add_argument("-m", "--minimumsize",
                    help="Minimum size (in arcsec) for the regions to be considered for the search",
                    type=float, required=False, default=15.0)

parser.add_argument("-s", "--significance",
                    help="Significance threshold for the search, in sigma",
                    type=float, required=True, default=5.0)

parser.add_argument("-x", "--expomap",
                    help="Exposure map for whole observation",
                    type=validFITSfile, required=True)

# parser.add_argument("-c","--catalog",
#                    help="FITS file containing the source catalog",
#                    type=validFITSfile,required=True)

# parser.add_argument("-o","--outfile",
#                    help="Name for the output file containing the list of found transients",
#                    required=True)

parser.add_argument("-v", "--verbosity", required=False, default='info',
                    choices=['info', 'debug'])


def go(args):
    # Set up the logger
    # Set up the logger
    levels = {'info': logging.INFO, 'debug': logging.DEBUG}
    log.setLevel(levels[args.verbosity])

    log.info("Probing the Hardware Unit...")
    hwu = HardwareUnit.hardwareUnitFactory(args.eventfile)

    # Root for the name of the files
    root = ".".join(os.path.basename(args.eventfile).split(".")[:-1])

    log.info("Analyzing %s" % (hwu.getName()))

    # Read in the regions

    log.info("Reading regions and filtering events...")

    with warnings.catch_warnings(record=True) as w:

        # Cause all warnings to be ignored

        warnings.simplefilter("ignore")

        ff = fitsRegions.FitsRegionFile(args.eventfile, args.minimumsize / hwu.getPixelScale())

    log.info("done")

    # Find the excesses in the regions

    log.info("Looking for excesses in regions...")

    allExcesses = []

    for i, (x, y, t, regfilt, region) in enumerate(ff.iteritems()):

        log.info("Processing region %s" % region)

        if len(x) < 3:

            # Empty region. Skip it
            log.info("Region has less than 3 counts. Skip it." % region)

            continue

        else:

            log.info("Region has %s events." % (len(x)))

        fbs = FixedBinSearch.FixedBinSearch(x, y, t,
                                            args.dt,
                                            args.eventfile,
                                            args.expomap,
                                            hwu,
                                            region)

        excesses, figs = fbs.scan(regfilt, args.significance)

        allExcesses.extend(excesses)

        for j, fig in enumerate(figs):
            fig.savefig("%s_r%s_%s.png" % (root, i + 1, j + 1), tight_layout=True)

            # Explicitly close the figure to save memory
            plt.close(fig)

        log.info("done")

    log.info("%s excesses found before consolidation" % len(allExcesses))

    allExcesses = numpy.asarray(allExcesses)

    log.debug("Excesses before consolidation:")

    for exc in allExcesses:
        log.debug("%.2f %.2f %.2f %.2f %.1f" % tuple(exc))

    finalExcesses = FixedBinSearch.consolidateTimeIntervals(allExcesses, hwu.getPixelScale())

    log.info("Found %s candidate transient episodes" % (len(finalExcesses)))

    # Get the WCS to transform from X,Y to RA,Dec
    wcs = XMMWCS.XMMWCS(args.eventfile)

    with open("%s_all.txt" % (root), "w+") as f:

        f.write("#RA Dec Tstart Tstop Probability X Y\n")

        for exc in finalExcesses:
            t1, t2, x, y, sigma = exc

            ra, dec = wcs.xy2sky([[x, y]])[0]

            # One-sided probability of a 5 sigma effect
            prob = scipy.stats.norm().sf(sigma) / 2.0

            f.write("%.5f %.5f %.2f %.2f %.3g %s %s\n" % (ra, dec, t1, t2, prob, x, y))

    # Second output: unique sources
    # Compute the distance between the reference
    # and all the others

    uniqueSources = []

    for exc in finalExcesses:

        t1, t2, x, y, sigma = exc

        # Compute the distance between this exc and all the others

        d = numpy.sqrt(numpy.power(x - finalExcesses[:, 2], 2) +
                       numpy.power(y - finalExcesses[:, 3], 2))

        # Note that this selection include "exc" itself (d is zero)
        idx = (d <= 10.0 / hwu.getPixelScale()) & (finalExcesses[:, -1] > 0)

        if numpy.sum(idx) > 0:
            # Among all the intervals selected, choose the one with the highest
            # significance

            iid = finalExcesses[idx, -1].argmax()

            uniqueSources.append(finalExcesses[idx][iid])

            # Now mark all these intervals by putting the sigma to -1

            finalExcesses[idx, -1] = -1.0

    log.info("Found %s unique new sources" % (len(uniqueSources)))

    with open("%s_unique.txt" % (root), "w+") as f:

        f.write("#RA Dec Tstart Tstop Probability\n")

        for exc in uniqueSources:
            t1, t2, x, y, sigma = exc

            ra, dec = wcs.xy2sky([[x, y]])[0]

            # One-sided probability of a n sigma effect
            prob = scipy.stats.norm().sf(sigma) / 2.0

            f.write("%.4f %.4f %.3f %.3f %.3g\n" % (ra, dec, t1, t2, prob))

    log.debug("all done")


# Main code
if __name__ == "__main__":
    # Parse arguments
    args = parser.parse_args()
    go(args)
