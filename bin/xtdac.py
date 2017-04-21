#!/usr/bin/env python

# Author: Giacomo Vianello (giacomov@stanford.edu)

##########################
# Plotting configuration
##########################

from XtDac.DivideAndConquer import matplotlibConfig

# Use a matplotlib backend which does not show plots to the user
# (they will be saved in files)
import matplotlib

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
import warnings
import math
import logging
import multiprocessing

import scipy.stats

import time as timemod

try:
    import astropy.io.fits as pyfits
except:
    # If this fail there is no way out
    import pyfits
pass

from XtDac.DivideAndConquer import GridGen
from XtDac.DivideAndConquer import HardwareUnit
from XtDac.DivideAndConquer import InterestingRegion
from XtDac.DivideAndConquer import Results
from XtDac.DivideAndConquer import TimeIntervalConsolidator
from XtDac.DivideAndConquer import XMMWCS
from XtDac.DivideAndConquer import Box

from XtDac.FixedBinSearch import Likelihood
from XtDac.FixedBinSearch import fitsRegions

time, X, Y, tstart, tstop, event_header = (None, None, None, None, None, None)

# Set up the logger (its verbosity level will be changed later on
# using the value provided by the user)
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("XtDac")


def validFITSfile(arg):

    if not os.path.isfile(arg):

        log.error("The file %s does not exist!" % arg)
        sys.exit()

    # Try to open it to see if it is a FITS file
    try:

        fits_file = pyfits.open(arg)

    except:

        log.error("The file %s is not a valid FITS file" % arg)
        sys.exit()

    else:

        fits_file.close()

    return arg


# def _initProcess(count, eventFile):
#
#     global counter
#     global time, X, Y, tstart, tstop
#
#     counter = count
#
#     # Read these here, so they will be read only one time for each of the
#     # processes spawn by multiprocessing
#     with pyfits.open(eventFile) as fitsFile:
#
#         time = fitsFile['EVENTS'].data.field("TIME")
#         X = fitsFile['EVENTS'].data.field("X")
#         Y = fitsFile['EVENTS'].data.field("Y")
#
#         # Sort by time
#         idx = numpy.argsort(time)
#         time = time[idx]
#         X = X[idx]
#         Y = Y[idx]
#
#         tstart = fitsFile['EVENTS'].header.get("TSTART")
#         tstop = fitsFile['EVENTS'].header.get("TSTOP")
#

def trueWorker(box_def, eventFile, nullHypProb, totRegions, bkgIdx=None):

    box = Box.Box(*box_def)

    box.setEventfile(eventFile)

    box.readEvents(True, time, X, Y, tstart, tstop)

    if (bkgIdx is not None):

        # Use custom background (usually when looking in regions close to
        # bright sources)

        box.setBackground(bkgIdx)


    # sys.stderr.write(str(box.nOutOfRegion))

    if (box.isEmpty() or box.nInRegion < 2 or box.filledArea < 0.5):

        results = []

    else:

        results = box.findExcesses(nullHypProb)

    # box.clearMemory()

    return results


def array2string(array):
    return ",".join(map(lambda x: "%s" % x, array))


# from memory_profiler import profile

def _process_regions(eventfile, typeIerror, regions, n_cpus, log):

    # Define a wrapper for the trueWorker to avoid duplicating all arguments
    workerWrapper = functools.partial(trueWorker,
                                      eventFile=eventfile,
                                      nullHypProb=typeIerror,
                                      totRegions=len(regions),
                                      bkgIdx=None)

    # Send each region to one worker to be processed,
    # and collect the results

    n_regions = len(regions)

    results = []

    # This is the unit for reporting progress

    chunk_size = 30

    progress_unit = max(chunk_size * n_cpus, n_regions / 20)

    #progress_unit = max(n_cpus, int(float(n_regions) / 10.0))

    if n_cpus > 1:

        # Parallel version

        log.debug("Starting a pool of %s python processes" % n_cpus)

        pool = multiprocessing.Pool(n_cpus,
                                    # initializer=_initProcess,
                                    # initargs=(counter, args.eventfile)
                                    )

        log.debug("Feeding the regions to the processes...")

        for i, result in enumerate(pool.imap(workerWrapper, regions, chunksize=chunk_size)):

            if (i + 1) % progress_unit == 0 or (i + 1) == n_regions:

                log.info("%s out of %s regions completed (%.0f percent)" % (i + 1, n_regions,
                                                                              (i + 1) / float(n_regions) * 100))

            results.append(result)

        pool.close()
        pool.join()

    else:

        # Serial version
        log.debug("Using serial version (only 1 CPU)")

        for i, region in enumerate(regions):

            results.append(workerWrapper(region))

            if (i + 1) % progress_unit == 0 or (i + 1) == n_regions:

                log.info("%s out of %s regions completed (%.0f percent)" % (i + 1, n_regions,
                                                                              (i + 1) / float(n_regions) * 100))

    assert len(results) == n_regions, "Something went wrong in the computation. The number of results doesn't match " \
                                      "the number of regions"

    return results

# @profile
def go(args):

    global time, X, Y, tstart, tstop, event_header

    # Set up the logger
    levels = {'info': logging.INFO, 'debug': logging.DEBUG}
    log.level = levels[args.verbosity]

    log.debug("Command line: %s" % (" ".join(sys.argv)))

    # Instantiate the HardwareUnit class
    log.debug("Probing Hardware Unit...")
    hwu = HardwareUnit.hardwareUnitFactory(args.eventfile)
    log.debug("probed %s" % hwu.getName())
    log.debug("done")

    log.debug("Reading event file...")

    # Remember:  X, Y, tstart, tstop, time, event_header are global variables (module-level)

    with pyfits.open(args.eventfile) as fitsFile:

        # Get events

        X_ = fitsFile['EVENTS'].data.field("X")
        Y_ = fitsFile['EVENTS'].data.field("Y")
        time_ = fitsFile['EVENTS'].data.field("TIME")

        # Order them by time

        # Note that this is advanced slicing, hence it returns a COPY of X_, Y_ and time_.
        # That's why we delete them explicitly immediately after, to save memory

        idx = numpy.argsort(time_)
        time = time_[idx]
        X = X_[idx]
        Y = Y_[idx]

        event_header = fitsFile['EVENTS'].header

        # Get the start of the first GTI and the stop of the last one
        gti_starts = []
        gti_stops = []
        for ext in fitsFile[1:]:

            if ext.header['EXTNAME'] == 'GTI':

                gti_starts.append(ext.data.field("START").min())
                gti_stops.append(ext.data.field("STOP").max())

        # Since we might have randomized the times (in the Chandra pipeline,  in the xtc_filter_event_file script),
        #  we need to make sure that the tstart is either the beginning of the first GTI or the time of the first event

        tstart = min(min(gti_starts), time_.min() - 1e-3)
        tstop = max(max(gti_stops), time_.max() + 1e-3)

        # Now make arrays read-only so they will never be copied

        X.setflags(write=False)
        Y.setflags(write=False)
        time.setflags(write=False)

        # Save memory

        del X_, Y_, time_


    log.debug("done")

    # Instantiate the GridGen class and generate the grid
    # of regions
    log.debug("Generating grid...")

    if hwu.getName().find("ACIS")==0:

        # Chandra. Automatically compute the step size, as the average
        # size of the PSF in the detector

        try:

            import psf
            import caldb4
            from astropy.coordinates import SkyCoord
            import astropy.units as u

        except ImportError:

            raise RuntimeError("Cannot import psf and/or caldb4 and or astropy module from CIAO. "
                               "Is CIAO python configured?")

        cdb = caldb4.Caldb(telescope="CHANDRA", product="REEF")
        reef = cdb.search[0]
        extno = cdb.extno()
        # Replace the trailing '[..]' block number specifier
        reef = reef.split('[')[0] + "[{}]".format(extno + 1)

        pdata = psf.psfInit(reef)

        # Get the coordinates of the events at the corners of the CCD
        # (don't use the minimum and maximum of X and Y because the CCD is rotated
        # with respect to sky coordinates)

        # Get the corners

        minx = X.min()
        maxx = X.max()
        miny = Y.min()
        maxy = Y.max()

        corner_1 = [minx, Y[X == minx].max()]
        corner_2 = [X[Y == maxy].max(), maxy]
        corner_3 = [X[Y == miny].max(), miny]
        corner_4 = [maxx, Y[X == maxx].max()]

        # Get the aim point
        ra_pointing = event_header.get('RA_PNT')
        dec_pointing = event_header.get('DEC_PNT')

        system = event_header.get("RADECSYS")

        if system is None:

            system = 'ICRS'

        wcs = XMMWCS.XMMWCS(args.eventfile, X, Y)

        x_pointing, y_pointing = wcs.sky2xy([[ra_pointing, dec_pointing]])[0]

        # Computing maximum distance between corners and the pointing
        distances_to_corners = numpy.linalg.norm(numpy.array([x_pointing, y_pointing]) -
                                                 numpy.array([corner_1, corner_2, corner_3, corner_4]),
                                                 axis=1)

        max_distance = distances_to_corners.max()

        pointing = SkyCoord(ra=ra_pointing * u.degree, dec=dec_pointing * u.degree, frame=system.lower())

        def get_theta(x_, y_):

            point = wcs.xy2sky([[x_, y_]])[0]

            c1 = SkyCoord(ra=point[0] * u.degree, dec=point[1] * u.degree, frame=system.lower())

            this_theta = c1.separation(pointing)

            # Return it in arcmin

            return this_theta.to(u.arcmin).value

        def get_psf_size(x_, y_):

            theta = get_theta(x_, y_)

            return psf.psfSize(pdata, 1.5, theta, 0.0, 0.68)

        gg = GridGen.GridGenChandra(hwu)

        gg.makeGrid(x_pointing, y_pointing, get_psf_size, max_distance)

    else:

        # Something else. Use provided step size

        gg = GridGen.GridGen(hwu)

        gg.makeGrid(args.regionsize, 'sky', args.multiplicity)

    log.debug("done")

    # Get the boxes definitions

    regions = gg.getBoxes()

    # with open('test_variable_size.reg','w+') as f:
    #
    #     for spec in regions:
    #
    #         reg = Box.Box(*spec)
    #
    #         f.write("%s" % "".join(reg.getDs9Region()))

    #sys.exit(0)

    n_events = time.shape[0]

    log.info("Processing interval %s - %s (duration: %s s, %s events)"
             % (tstart, tstop, tstop - tstart, n_events))

    results = _process_regions(args.eventfile, args.typeIerror, regions, args.ncpus, log)

    log.debug("done")

    log.debug("Selecting interesting regions with more than one interval...")

    interesting_regions = []

    for i, intervals in enumerate(results):

        # Remember: intervals contains the edges of the intervals

        if len(intervals) > 2:

            box = Box.Box(*regions[i])

            box.setEventfile(args.eventfile)

            box.readEvents(preLoaded=True, time=time, X=X, Y=Y, tstart=tstart, tstop=tstop)

            #if args.writeRegionFiles == 'yes':

            #    box.writeRegion("%s_exc%i.reg" % (root_filename, i + 1))
            #    box.writeRegion("%s_exc%i.reg" % (root_filename, i + 1))

            log.debug("Region %i has %i intervals" % (i+1, len(intervals)-1))

            box.writeRegion("interesting_region_%i.reg" % (i+1))

            for j,(t1,t2) in enumerate(zip(intervals[:-1],intervals[1:])):

                log.debug("  %i : %s - %s (%s s long)" % (j+1, t1, t2, t2-t1))

            thisInterestingRegion = InterestingRegion.InterestingRegion(box, intervals,
                                                                        time, X, Y, tstart, tstop)
            interesting_regions.append(thisInterestingRegion)

    log.debug("done")

    log.debug("Removing overlapping intervals...")

    consolidator = TimeIntervalConsolidator.TimeIntervalConsolidator(interesting_regions)

    cleanedIntervals = consolidator.consolidate()

    log.debug("done")

    log.info("Kept %s intervals" % len(cleanedIntervals))

    if args.sigmaThreshold > 0:

        # Now for each cleaned interval perform a likelihood analysis on a region
        # larger than the box. This is needed otherwise it is too difficult to distinguish
        # a PSF-like excess and a flat-like excess

        log.info("Computing final significance...")

        finalCandidates = []

        for i, interval in enumerate(cleanedIntervals):

            log.info("Processing interval %s of %s" % (i + 1, len(cleanedIntervals)))

            this_interval_duration = interval.tstop - interval.tstart

            log.info("duration: %.1f s" % this_interval_duration)

            if this_interval_duration > args.max_duration:

                log.info("Duration is longer than the maximum allowed of %.1f" % (args.max_duration))

                continue

            if interval.nEvents < args.min_number_of_events:

                log.info("Less than %s events, discarding candidate with %s events" % (args.min_number_of_events,
                                                                                       interval.nEvents))

                continue

            boxx, boxy = interval.interestingRegion.getCenterPhysicalCoordinates()

            boxw, boxh = interval.interestingRegion.box.width, interval.interestingRegion.box.height

            # Search region is 25 arcsec of radius

            if hwu.getName().find("ACIS")==0:

                # For Chandra, let's adapt to the size of the PSF (but keep a minimum size of
                # at least 50 px)

                box_size_arcsec = max([boxw, boxh]) * hwu.getPixelScale()

                radius = max(50.0, (box_size_arcsec /2.0) * 1.5 / hwu.getPixelScale())

                log.info("Radius for search region: %s px" % radius)

            else:

                radius = 40.0 / hwu.getPixelScale()

            searchRegionStr = 'physical;circle(%s,%s,%s)' % (boxx, boxy, radius)

            idx = (time >= interval.tstart) & (time <= interval.tstop)

            assert X[idx].shape[0] > 0
            assert Y[idx].shape[0] > 0

            with warnings.catch_warnings():

                # Cause all warnings to be ignored

                warnings.simplefilter("ignore")

                searchRegionDef = fitsRegions.FitsRegion(X[idx], Y[idx], time[idx],
                                                         event_header, searchRegionStr)

            # This is actually a loop of only one iteration

            for x, y, t, regfilter, reg in searchRegionDef.iteritems():

                xmin = boxx - radius
                xmax = boxx + radius
                ymin = boxy - radius
                ymax = boxy + radius

                if hwu.getName().find("ACIS") == 0:

                    # For Chandra, let's adapt to the size of the PSF, but make at least
                    # twice the size of the pixel (remember, this is in pixels)

                    r_binsize = max(radius / 60.0, 2.0)

                else:

                    r_binsize = 40.0 / hwu.getPixelScale() / 10.0

                assert (xmax-xmin) / r_binsize > 3.0, "Not enough bins for likelihood!"

                searchRegion = Likelihood.Region(xmin, xmax,
                                                 ymin, ymax,
                                                 r_binsize, regfilter,
                                                 args.expomap, args.eventfile)

                ls = Likelihood.Likelihood(x, y, searchRegion)

            # Build the likelihood model

            # Bkg model (isotropic)
            iso = Likelihood.Isotropic("bkg", 1.0)

            m = Likelihood.GlobalModel("likeModel") + iso

            try:

                ls.setModel(m)

            except UnboundLocalError:

                import pdb;pdb.set_trace()

            # Minimize the mlogLike to get the background level
            like0 = ls.minimize(verbose=0)

            # Small TS map to figure out the source position

            trial_x = numpy.linspace(boxx - boxw / 2.0, boxx + boxw / 2.0, 16)
            trial_y = numpy.linspace(boxy - boxh / 2.0, boxy + boxh / 2.0, 16)

            ra, dec = interval.interestingRegion.getCenterSkyCoordinates()

            maxTS = 0
            _best_x = None
            _best_y = None

            TSmap = numpy.zeros([trial_y.shape[0], trial_x.shape[0]])

            # import cProfile, pstats, StringIO
            # pr = cProfile.Profile()
            # pr.enable()

            for k, srcx in enumerate(trial_x):

                for h, srcy in enumerate(trial_y):

                    sys.stderr.write(".")

                    # Now fit adding a source at the position of the maximum of the
                    # image of this region

                    # srcx, srcy = ls.getMaximumPosition()

                    if h == 0:

                        # This is the very first loop. Make the image of the PSF.
                        # We assume that the PSF does not change much within the region, so we
                        # will use the same image for all the next iterations

                        thisSrc = Likelihood.PointSource("testSrc", srcx, srcy, args.eventfile, hwu)
                        psffile = thisSrc.outfile

                    else:

                        # Re-use the psf file already computed
                        thisSrc = Likelihood.PointSource("testSrc", srcx, srcy,
                                                         args.eventfile, hwu,
                                                         psffile)

                    pass

                    m += thisSrc

                    ls.setModel(m)

                    like1 = ls.minimize(verbose=0)

                    TS = 2 * (like0 - like1)

                    TSmap[h, k] = TS

                    if TS > maxTS:

                        maxTS = TS

                        _best_x = srcx
                        _best_y = srcy

                        with warnings.catch_warnings():

                            # Cause all warnings to be ignored

                            warnings.simplefilter("ignore")

                            ls.plot().savefig("c_%.4f_%.4f_%i.png" % (ra,dec,i))

                    m.removeSource("testSrc")

            # pr.disable()
            # s = StringIO.StringIO()
            # sortby = 'cumulative'
            # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            # ps.print_stats()
            # print s.getvalue()
            # sys.exit(0)

            sys.stderr.write("\n")

            significance = math.sqrt(max(TSmap.max(), 0))

            log.info("number of events: %s" % interval.nEvents)

            log.info("Region @ (RA,Dec) = (%.4f,%.4f) (%.3f - %.3f s) "
                     "-> %.1f sigma" % (ra, dec,
                                        interval.tstart,
                                        interval.tstop,
                                        significance))

            if significance >= args.sigmaThreshold:

                log.debug("Keeping candidate")

                # One-sided probability of a n sigma effect
                prob = scipy.stats.norm().sf(significance) / 2.0

                interval.probability = prob

                interval._best_localization = (_best_x, _best_y)

                finalCandidates.append(interval)

            else:

                log.debug("Discarding candidate")

                # plt.imshow(TSmap, interpolation='none', origin='lower')
                # plt.savefig("exc%s.png" % i, tight_layout=True)



        log.debug("done")

    else:

        log.debug("Threshold for final significance is <= 0, no likelihood analysis will be performed")

        finalCandidates = cleanedIntervals

        # The probability has not been computed, so let's assign -1
        for c in finalCandidates:
            c.probability = -1

    pass

    log.debug("Preparing the summary Ascii file...")

    root_filename = ".".join(os.path.basename(args.eventfile).split(".")[:-1])

    # I need to read the WCS from the event file, if transient_pos is true
    wcs = None

    if args.transient_pos:

        wcs = XMMWCS.XMMWCS(args.eventfile, X, Y)

    with open("%s_res.txt" % root_filename, "w+") as f:

        f.write("#RA Dec Tstart Tstop Probability\n")

        for i, interval in enumerate(finalCandidates):

            if args.transient_pos:

                # Write the position of the transient
                points = wcs.xy2sky([list(interval._best_localization)])

                ra, dec = points[0]

            else:

                # Write the position of the center of the region containing the transient

                ra, dec = interval.interestingRegion.getCenterSkyCoordinates()


            # Find the first and last event in the interval
            idx = (time >= interval.tstart) & (time <= interval.tstop)

            first_event_timestamp = time[idx].min()
            last_event_timestamp = time[idx].max()

            f.write("%.4f %.4f %.3f %.3f %.3g\n" % (ra, dec, first_event_timestamp - 1e-2, last_event_timestamp + 1e-2,
                                                    interval.probability))

            if args.writeRegionFiles == 'yes':

                interval.interestingRegion.box.writeRegion("%s_candidate_%i.reg" % (root_filename, i + 1))

    log.debug("done")

    jobStop = timemod.time()

    if args.html_summary:

        log.debug("Preparing the summary HTML file...")
        # Summarize the results in a HTML file
        resultsSummary = Results.Summary()
        resultsSummary.addJobInfo(jobStop - jobStart, args.ncpus, args.typeIerror)
        resultsSummary.addObsInfo(args.eventfile, hwu.getName(), tstart, tstop, n_events)
        resultsSummary.addWholeHWUlightCurve(time, tstart, tstop, figsize=(8, 8))
        resultsSummary.addResults(map(lambda interval: interval.interestingRegion, finalCandidates))
        root_filename = ".".join(os.path.basename(args.eventfile).split(".")[:-1])
        resultsSummary.write("%s_res.html" % root_filename)
        log.debug("done")

pass

if __name__ == "__main__":

    jobStart = timemod.time()

    # Define and parse input arguments

    desc = '''EXTraS Divide-and-Conquer algorithm to find
                                     transients.'''
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-e", "--eventfile",
                        help="Event file to use (already cleaned and selected in" +
                             " energy and quadrant/CCD)",
                        required=True, type=validFITSfile)

    parser.add_argument("-x", "--expomap",
                        help="Exposure map for the whole observation. It is only" +
                             " needed to figure out gaps and masked pixels",
                        required=True, type=validFITSfile)

    parser.add_argument("-b", "--backgroundRegion", help="Custom background circle. To be" +
                                                         " specified as ra,dec,radius. For example '-b 187.5," +
                                                         "-23.2,10' corresponds to a circle with 10 arcsec " +
                                                         "radius centered on R.A., Dec. = (187.5, -23.2)",
                        required=False, default=None)

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

    parser.add_argument("-s", "--sigmaThreshold",
                        help="Threshold for the final significance. All intervals found " +
                             "by the bayesian blocks " +
                             "algorithm which does not surpass this threshold will not be saved in the " +
                             "final file.",
                        type=float,
                        default=5.0,
                        required=False)

    parser.add_argument("-r", "--regionsize",
                        help="Approximate side for the square regions in which the" +
                             " search will be performed (Default: 60 arcsec)",
                        type=float,
                        default=60,
                        required=False)

    parser.add_argument("--max_duration",
                        help="Do not consider candidate transients with a duration longer then max_duration",
                        type=float,
                        default=1e9,
                        required=False)

    parser.add_argument("--min_number_of_events",
                        help="Do not consider candidate transients with less than this number of events",
                        type=int,
                        default=0,
                        required=False)


    parser.add_argument("-w", "--writeRegionFiles",
                        help="Write a ds9 region file for each region with excesses?",
                        type=str,
                        default='yes',
                        required=False,
                        choices=['yes', 'no'])

    parser.add_argument('--transient_pos', dest='transient_pos', action='store_true')
    parser.set_defaults(transient_pos=False)

    parser.add_argument('--no-html', dest='html_summary', action='store_false')
    parser.set_defaults(html_summary=True)

    parser.add_argument("-v", "--verbosity",
                        required=False,
                        default='info',
                        choices=['info', 'debug'])

    args = parser.parse_args()

    go(args)