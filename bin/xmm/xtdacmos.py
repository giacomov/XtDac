#!/usr/bin/env python

import argparse
import subprocess

import os

try:
    import astropy.io.fits as pyfits
except:
    # If this fail there is no way out
    import pyfits
pass
import logging
import sys

# Set up the logger (its verbosity level will be changed later on
# using the value provided by the user)
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger("xtdacmos")


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


def filterFITS(fitsfile, expr, output):
    if (os.path.exists(output)):
        os.remove(output)

    cmdLine = ("""fcopy '%s[EVENTS][%s]' %s"""
               % (fitsfile, expr, output))

    log.debug(cmdLine)
    subprocess.call(cmdLine, shell=True)


def callxtdac(filename, args):
    path = os.path.dirname(os.path.abspath(__file__))

    xtdacpath = os.path.join(path, "XtDac.py")

    cmdLine = ('%s -e %s -c %s -p %s -r %s -w %s -v %s'
               % (xtdacpath, filename, args.ncpus,
                  args.typeIerror, args.regionsize,
                  args.writeRegionFiles, args.verbosity))
    log.debug(cmdLine)
    subprocess.call(cmdLine, shell=True)


if __name__ == "__main__":

    # Define and parse input arguments

    desc = '''Apply EXTraS Divide-and-Conquer algorithm
                                     to a whole MOS dataset.'''
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-e", "--eventfile",
                        help="Event file to use (already cleaned)",
                        required=True, type=validFITSfile)

    parser.add_argument("-c", "--ncpus", help="Number of CPUs to use (default=1)",
                        type=int, default=1, required=False)

    parser.add_argument("-p", "--typeIerror",
                        help="Type I error probability for the Bayesian Blocks " +
                             "algorithm.",
                        type=float,
                        default=1e-6,
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
                        default='yes',
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

    bandsPI = [200, 500, 1000, 2000, 4500, 12000]
    ccds = range(1, 8)
    rootFile = ".".join(os.path.basename(args.eventfile).split(".")[:-1])

    # Divide the event files and run XtDac
    for ccd in ccds:

        log.info("CCD %i:" % (ccd))

        thisQFilename = "%s_ccd%i.fits" % (rootFile, ccd)

        filterFITS(args.eventfile,
                   "CCDNR == %s " % (ccd),
                   thisQFilename)

        # Count the number of events
        with pyfits.open(thisQFilename) as f:

            if (f['EVENTS'].data.shape[0] == 0):
                log.warning("\nNo event in this CCD!\n")
                continue

        log.info("Whole event file...")
        callxtdac(thisQFilename, args)

        for i, (emin, emax) in enumerate(zip(bandsPI[:-1], bandsPI[1:])):
            log.info("Band %i..." % (i + 1))

            thisFilename = "%s_ccd%i_b%i.fits" % (rootFile, ccd, i + 1)

            filterFITS(thisQFilename,
                       "PI >= %s && PI < %s" % (emin, emax),
                       thisFilename)

            callxtdac(thisFilename, args)
            log.info("done")
        pass
    pass
