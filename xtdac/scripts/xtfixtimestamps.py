#!/usr/bin/env python

# Author: Giacomo Vianello (giacomov@stanford.edu)

##############
# Pre-processing to account for duplicated time stamps in event files
##############

# Set up the logger (its verbosity level will be changed later on
# using the value provided by the user)
import logging

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)-15s %(message)s")

log = logging.getLogger("main")

import argparse
import os, sys
import numpy
import subprocess

try:

    import astropy.io.fits as pyfits

except:

    import pyfits

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


parser = argparse.ArgumentParser(description="Pre-processing of the event file to get rid of duplicated time stamps")

parser.add_argument("-e", "--eventfile",
                    help="Event file to use (already cleaned and selected in" +
                         " energy and quadrant/CCD), and containing the region definitions",
                    required=True, type=validFITSfile)

parser.add_argument("-v", "--verbosity", required=False, default='info',
                    choices=['info', 'debug'])


def go(args):

    # Set up the logger

    levels = {'info': logging.INFO, 'debug': logging.DEBUG}
    log.setLevel(levels[args.verbosity])

    log.debug("Opening FITS file")

    with pyfits.open(args.eventfile) as f:

        time = f["EVENTS"].data.TIME[:]

        idx = time.argsort()

        if not numpy.all(time[idx] == f["EVENTS"].data.field("TIME")):

            log.info("Events are out of order! Ordering them...")

            cmd_line = "fsort %s[EVENTS] TIME heap" % args.eventfile

            log.info("Executing:")
            log.info(cmd_line)

            subprocess.check_call(cmd_line, shell=True)

    u, uidx = numpy.unique(time, return_index=True)

    if u.shape[0] != time.shape[0]:

        log.info("Found %s events with the same time stamp. Correcting..." % (2*(time.shape[0] - u.shape[0])))

        # Re-open file and fix duplicated events

        with pyfits.open(args.eventfile, mode='update') as f:

            idx = numpy.where(time[1:] == time[:-1])[0]

            deltas = numpy.random.uniform(-1e-3,1e-3,idx.shape[0])

            deltas.sort()

            time[idx] += deltas

            f["EVENTS"].data.TIME[:] = time

    else:

        log.info("No events have the same time stamp.")


# Main code
if __name__ == "__main__":
    # Parse arguments
    args = parser.parse_args()
    go(args)