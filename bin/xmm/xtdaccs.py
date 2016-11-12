#!/usr/bin/env python

# Author: Giacomo Vianello (giacomov@stanford.edu)

##########################
# Plotting configuration
##########################

import matplotlib

from XtDac.bin.xtdac import matplotlibConfig

matplotlib.use("Agg")

rcParams = matplotlibConfig.getConfig()
for k, v in rcParams.iteritems():
    matplotlib.rcParams[k] = v
###########################

import argparse
import os
import sys
import functools
import numpy
import time as timemod

try:
    import astropy.io.fits as pyfits
except:
    # If this fail there is no way out
    import pyfits
pass

import logging
import multiprocessing

from XtDac.bin.xtdac import GridGen
from XtDac.bin.xtdac import HardwareUnit
from XtDac.bin.xtdac import InterestingRegion
from XtDac.bin.xtdac import Results
from XtDac.bin.xtdac import TimeIntervalConsolidator
from XtDac.bin.xtdac import XMMWCS

# These global variables will be used by the processes spawn by
# the multiprocessing module
counter = multiprocessing.Value('i', 0)
time, X, Y, tstart, tstop = (None, None, None, None, None)

# Set up the logger (its verbosity level will be changed later on
# using the value provided by the user)
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("xtdaccs")


def validFITSfile(arg):
    if not os.path.isfile(arg):
        log.error("The file %s does not exist!" % arg)
        sys.exit()

    # Try to open it to see if it is a FITS file
    try:
        pyfits.open(arg)
    except:
        log.error("The file %s is not a valid FITS file" % arg)
        sys.exit()
    return arg


def _initProcess(count, eventFile):
    global counter
    global time, X, Y, tstart, tstop

    counter = count

    # Read these here, so they will be read only one time for each of the
    # processes spawn by multiprocessing
    with pyfits.open(eventFile) as fitsFile:
        time = fitsFile['EVENTS'].data.field("TIME")
        X = fitsFile['EVENTS'].data.field("X")
        Y = fitsFile['EVENTS'].data.field("Y")

        # Sort by time
        idx = numpy.argsort(time)
        time = time[idx]
        X = X[idx]
        Y = Y[idx]

        tstart = fitsFile['EVENTS'].header.get("TSTART")
        tstop = fitsFile['EVENTS'].header.get("TSTOP")


def trueWorker(box, eventFile, nullHypProb, totRegions, bkgIdx=None):
    # sys.stderr.write(">")

    box.setEventfile(eventFile)

    box.readEvents(True, time, X, Y, tstart, tstop)

    if (bkgIdx is not None):
        # Use custom background (usually when looking in regions close to 
        # bright sources)

        box.setBackground(bkgIdx)

    pass

    # sys.stderr.write(str(box.nOutOfRegion))

    if (box.isEmpty() or box.nInRegion < 2 or box.filledArea < 0.5):
        # print("Region %s is empty" %(counter.value))
        results = []

    else:
        results = box.findExcesses(nullHypProb)

    with counter.get_lock():
        counter.value += 1

    box.clearMemory()

    sys.stderr.write("\r%s out of %s regions completed (%.0f percent)" % (counter.value, totRegions,
                                                                          counter.value / float(totRegions) * 100))

    return [box, results]


pass


def array2string(array):
    return ",".join(map(lambda x: "%s" % x, array))


if __name__ == "__main__":
    jobStart = timemod.time()

    # Define and parse input arguments

    desc = '''EXTraS Divide-and-Conquer algorithm to find
                                     transients CLOSE TO POINT SOURCES.'''
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-e", "--eventfile",
                        help="Event file to use (already cleaned and selected in" +
                             " energy and quadrant/CCD)",
                        required=True, type=validFITSfile)

    parser.add_argument("-s", "--sourcesfile", help="A FITS file containing RA,Dec of known sources",
                        required=True, type=validFITSfile)

    parser.add_argument("-m", "--multiplicity", help="Control the overlap of the regions." +
                                                     " A multiplicity of 2 means the centers of the regions are" +
                                                     " shifted by 1/2 of the region size (they overlap by 50 percent)," +
                                                     " a multiplicity of 4 means they are shifted by 1/4 of " +
                                                     " their size (they overlap by 75 percent), and so on.",
                        required=False, default=2.0, type=float)

    parser.add_argument("-c", "--ncpus", help="Number of CPUs to use (default=1)",
                        type=int, default=1, required=False)

    parser.add_argument("-p", "--typeIerror",
                        help="Type I error probability for the Bayesian Blocks " +
                             "algorithm.",
                        type=float,
                        default=1e-5,
                        required=False)

    parser.add_argument("-t", "--finalProbThreshold",
                        help="Threshold for the final probability. All intervals found " +
                             "by the bayesian blocks " +
                             "algorithm which does not surpass this threshold will not be saved in the " +
                             "final file.",
                        type=float,
                        default=1e-5,
                        required=False)

    parser.add_argument("-r", "--regionsize",
                        help="Approximate side for the square regions in which the" +
                             " search will be performed (Default: 60 arcsec)",
                        type=float,
                        default=60,
                        required=False)

    parser.add_argument("-w", "--writeRegionFiles",
                        help="Write a ds9 region file for each region with excesses?",
                        type=str,
                        default='no',
                        required=False,
                        choices=['yes', 'no'])

    parser.add_argument("-v", "--verbosity",
                        required=False,
                        default='info',
                        choices=['info', 'debug'])

    args = parser.parse_args()

    # Set up the logger
    levels = {'info': logging.INFO, 'debug': logging.DEBUG}
    log.level = levels[args.verbosity]

    # Instantiate the HardwareUnit class
    log.debug("Probing Hardware Unit...")
    hwu = HardwareUnit.hardwareUnitFactory(args.eventfile)
    log.debug("done")

    log.debug("Reading known sources...")

    with pyfits.open(args.sourcesfile) as f:

        known_sources_ = zip(f['SRCLIST'].data.field("RA"),
                             f['SRCLIST'].data.field("DEC"))

    log.debug("done")

    log.debug("Transforming RA,Dec in X,Y for known sources...")

    wcs = XMMWCS.XMMWCS(args.eventfile)
    known_sources = wcs.sky2xy(known_sources_)

    log.debug("done")

    # Instantiate the GridGen class and generate the grid
    # of regions
    log.debug("Generating grid...")
    gg = GridGen.GridGen(hwu)
    gg.makeGrid(args.regionsize, 'sky', args.multiplicity, known_sources)
    log.debug("done")

    regions = gg.getBoxes()

    # Open the event file and print some info
    log.debug("Reading events file...")
    with pyfits.open(args.eventfile) as fitsFile:
        time = fitsFile['EVENTS'].data.field("TIME")
        X = fitsFile['EVENTS'].data.field("X")
        Y = fitsFile['EVENTS'].data.field("Y")
        tstart = fitsFile['EVENTS'].header.get("TSTART")
        tstop = fitsFile['EVENTS'].header.get("TSTOP")
        Nevents = time.shape[0]

    log.debug("done")

    log.debug("Assigning closest source to each region...")
    # Assign the closest source to each region
    for box in regions:
        xc, yc = (box.c1 + box.width / 2.0, box.c2 + box.height / 2.0)

        dist2 = numpy.array(map(lambda (x, y): (xc - x) ** 2 + (yc - y) ** 2, known_sources))
        idx = dist2.argmin()

        box.setBrightSource(known_sources[idx][0], known_sources[idx][1])

    log.debug("done")

    log.info("Processing interval %s - %s (duration: %s s, %s events)"
             % (tstart, tstop, tstop - tstart, Nevents))

    # Define a wrapper for the trueWorker to avoid duplicating all arguments
    workerWrapper = functools.partial(trueWorker,
                                      eventFile=args.eventfile,
                                      nullHypProb=args.typeIerror,
                                      totRegions=len(regions),
                                      bkgIdx=None)

    # Send each region to one worker to be processed,
    # and collect the results
    log.debug("Starting a pool of %s python processes" % args.ncpus)
    pool = multiprocessing.Pool(args.ncpus,
                                initializer=_initProcess,
                                initargs=(counter, args.eventfile))
    # _initProcess(counter,eventFile)
    # results                       = map(workerWrapper, regions)
    log.debug("Feeding the regions to the processes...")
    results = pool.map(workerWrapper, regions)
    # This is to end the line in which I was reporting the progress
    sys.stderr.write("\n")

    log.debug("done")
    pool.close()

    log.debug("Selecting interesting regions with more than one interval...")
    interestingRegions = []
    rootFilename = ".".join(os.path.basename(args.eventfile).split(".")[:-1])

    for i, (box, intervals) in enumerate(results):

        if (len(intervals) > 2):

            # Sometimes at the beginning of the interval there are very short excesses
            # due to the integral distribution having large statistical uncertainties
            # Filter them out
            t1, t2 = intervals[:2]
            if (t2 - t1 <= 10):
                # Too short, filter it out
                cleanedIntervals = intervals[2:]
                if (len(cleanedIntervals) <= 2):
                    continue
                pass
            else:
                cleanedIntervals = intervals
            pass

            thisInterestingRegion = InterestingRegion.InterestingRegion(box, intervals,
                                                                        time, X, Y, tstart, tstop)
            interestingRegions.append(thisInterestingRegion)
        pass
    pass
    log.debug("done")

    log.debug("Removing overlapping intervals, and keep only excesses above floor probability...")
    consolidator = TimeIntervalConsolidator.TimeIntervalConsolidator(interestingRegions)
    cleanedIntervals = consolidator.consolidate(args.finalProbThreshold)
    log.debug("done")

    log.debug("Preparing the summary Ascii file...")

    rootFilename = ".".join(os.path.basename(args.eventfile).split(".")[:-1])

    cleanedRegions = []

    with open("%s_res.txt" % rootFilename, "w+") as f:

        f.write("#R.A. Dec Tstart Tstop Probability\n")

        for i, interval in enumerate(cleanedIntervals):
            ra, dec = interval.interestingRegion.getCenterSkyCoordinates()
            f.write("%.4f %.4f %.3f %.3f %.3g\n" % (ra, dec, interval.tstart, interval.tstop, interval.probability))

            cleanedRegions.append(interval.interestingRegion)

            if (args.writeRegionFiles == 'yes'):
                interval.interestingRegion.box.writeRegion("%s_exc%i.reg" % (rootFilename, i + 1))

        pass
    pass

    log.debug("done")

    jobStop = timemod.time()

    log.debug("Preparing the summary HTML file...")
    # Summarize the results in a HTML file
    resultsSummary = Results.Summary()
    resultsSummary.addJobInfo(jobStop - jobStart, args.ncpus, args.typeIerror)
    resultsSummary.addObsInfo(args.eventfile, hwu.getName(), tstart, tstop, Nevents)
    resultsSummary.addWholeHWUlightCurve(time, tstart, tstop, figsize=(8, 8))
    resultsSummary.addResults(cleanedRegions)
    rootFilename = ".".join(os.path.basename(args.eventfile).split(".")[:-1])
    resultsSummary.write("%s_res.html" % rootFilename)
    log.debug("done")

pass
