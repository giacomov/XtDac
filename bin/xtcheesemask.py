#!/usr/bin/env python

# Author: Giacomo Vianello (giacomov@stanford.edu)

##############
# Oversample an image and filter it according to a region file.
# Optionally, the output image can be a mask of 0 or 1
##############

import argparse
import logging

import astropy.io.fits as pyfits
import numpy
import os
import pyregion
import warnings
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.wcs.wcs import FITSFixedWarning
from scipy import ndimage

from XtDac.bin.xtdac import HardwareUnit
from XtDac.bin.xtdac import XMMWCS

levels = {'info': logging.INFO, 'debug': logging.DEBUG}
logging.basicConfig(level=levels['info'])
log = logging.getLogger("xtcheesemask")


class NotFITS(RuntimeError):
    pass


def resampleImage(imageData, imageWCS, scaleFactor, threshold):
    if type(scaleFactor) == int or type(scaleFactor) == float:
        scaleFactor = [float(scaleFactor), float(scaleFactor)]

        # Resample with constant interpolation
    mask = ndimage.zoom(imageData, scaleFactor, order=0, mode='nearest')

    # Make a mask
    idx = mask <= threshold
    mask[idx] = 0
    mask[~idx] = 1

    # Resample with linear interpolation
    scaledData = ndimage.zoom(imageData, scaleFactor, order=1, mode='nearest')

    # Restore zeros
    scaledData *= mask

    del mask

    # Take care of offset due to rounding in scaling image to integer pixel dimensions
    properDimensions = numpy.array(imageData.shape) * scaleFactor
    offset = properDimensions - numpy.array(scaledData.shape)

    # Rescale WCS
    try:
        oldCRPIX1 = imageWCS['CRPIX1']
        oldCRPIX2 = imageWCS['CRPIX2']
        CD11 = imageWCS['CD1_1']
        CD21 = imageWCS['CD2_1']
        CD12 = imageWCS['CD1_2']
        CD22 = imageWCS['CD2_2']
    except KeyError:
        # Try the older FITS header format
        try:
            oldCRPIX1 = imageWCS['CRPIX1']
            oldCRPIX2 = imageWCS['CRPIX2']
            CD11 = imageWCS['CDELT1']
            CD21 = 0
            CD12 = 0
            CD22 = imageWCS['CDELT2']
        except KeyError:
            scaledWCS = imageWCS.copy()
            return {'data': scaledData, 'wcs': scaledWCS}

    CDMatrix = numpy.array([[CD11, CD12], [CD21, CD22]], dtype=numpy.float64)
    scaleFactorMatrix = numpy.array([[1.0 / scaleFactor[0], 0], [0, 1.0 / scaleFactor[1]]])
    scaledCDMatrix = numpy.dot(scaleFactorMatrix, CDMatrix)

    scaledWCS = imageWCS.copy()
    scaledWCS['NAXIS1'] = scaledData.shape[1]
    scaledWCS['NAXIS2'] = scaledData.shape[0]
    scaledWCS['CRPIX1'] = oldCRPIX1 * scaleFactor[0] + offset[1]
    scaledWCS['CRPIX2'] = oldCRPIX2 * scaleFactor[1] + offset[0]
    scaledWCS['CD1_1'] = scaledCDMatrix[0][0]
    scaledWCS['CD2_1'] = scaledCDMatrix[1][0]
    scaledWCS['CD1_2'] = scaledCDMatrix[0][1]
    scaledWCS['CD2_2'] = scaledCDMatrix[1][1]

    return scaledData, scaledWCS


class Regions(object):
    def __init__(self, regionFile, headerWCS, reverse=False):

        self.reverse = bool(reverse)

        # Check if the file is a FITS file

        try:

            self._readFITSFile(regionFile, headerWCS)

        except NotFITS:

            # NO. Assume it is a ds9 region file

            self._readDs9File(regionFile, headerWCS)

    def _readFITSFile(self, regionFile, headerWCS):

        try:

            with pyfits.open(regionFile) as f:

                if 'REGION' in f:

                    # XMM-Newton file

                    data = f['REGION'].data

                elif 'SRCREG' in f:

                    # Chandra file

                    data = f['SRCREG'].data

        except:

            raise NotFITS("This is not a FITS file")

        else:

            hwu = HardwareUnit.hardwareUnitFactory(regionFile)
            xmmwcs = XMMWCS.XMMWCS(regionFile)

            lines = ['fk5\n']

            for i, row in enumerate(data):

                try:

                    X = row.field("X")[0]
                    Y = row.field("Y")[0]

                except IndexError:

                    # Probably a Chandra file
                    X = float(row.field("X"))
                    Y = float(row.field("Y"))

                ra, dec = xmmwcs.xy2sky([[X, Y]])[0]

                r1, r2 = row.field("R")[:2] * hwu.getPixelScale()  # arcsec
                rot = row.field("ROTANG")[0]

                reg = 'ellipse({ra},{dec},{r1}",{r2}",{rot})\n'.format(**locals())

                if reg.find("nan") >= 0:
                    # Sometimes some regions are invalid. Just ignore the problem
                    continue

                lines.append(reg)

            self._parseLines(lines, headerWCS)

    def _readDs9File(self, regionFile, headerWCS):

        with open(regionFile, "r") as f:
            lines = f.readlines()

        # Filter out comments, i.e., lines starting with '#'
        lines = filter(lambda x: x[0] != "#", lines)

        # print lines

        self._parseLines(lines, headerWCS)

    def _parseLines(self, lines, headerWCS):

        regDefinitions = filter(lambda x: x != 'fk5\n', lines)

        if len(regDefinitions) > 0:
            # Extract all circle(ra,dec,radius) expressions
            # (catch the warnings generated by pyregion, which are
            # harmless)

            with warnings.catch_warnings():

                warnings.filterwarnings(action="ignore", category=AstropyDeprecationWarning)
                warnings.filterwarnings(action="ignore", category=FITSFixedWarning)

                regions = pyregion.parse(" ".join(lines)).as_imagecoord(headerWCS)

            self.filter = regions

            self.n = len(regions)

        else:

            # There are no regions!

            self.filter = None

            self.n = 0

    def getMask(self, shape):

        msk = self.filter.get_mask(shape=shape)

        if self.reverse:
            # Reverse the selection (i.e., include regions instead of exclude them)

            msk = ~msk

        return msk

    def getNregions(self):
        return self.n


def validFITSfile(arg):
    existingFile(arg)

    # Try to open it to see if it is a FITS file
    try:
        pyfits.open(arg)
    except:
        msg = "The file %s is not a valid FITS file" % arg
        log.error(msg)
        raise argparse.ArgumentTypeError(msg)
    return arg


def existingFile(arg):
    if not os.path.isfile(arg):
        msg = "The file %s does not exist!" % arg
        log.error(msg)
        raise argparse.ArgumentTypeError(msg)
    return arg


def notExistingFile(arg):
    if os.path.isfile(arg):
        msg = "The output file %s already exists!" % arg
        log.error(msg)
        raise argparse.ArgumentTypeError(msg)
    return arg


parser = argparse.ArgumentParser(description="Filter an image according to a region file from ds9")

parser.add_argument("-i", "--image", help="Input image", type=validFITSfile, required=True)
parser.add_argument("-r", "--regionfile",
                    help="File containing the region definition (in ds9 format) or eventfile with REGION extension",
                    type=existingFile, required=True)
parser.add_argument("-m", "--makemask", help="If specified, the output will be a mask with only 0 and 1",
                    dest='makemask', action='store_true')
parser.add_argument("-o", "--outfile", help="Name for the output image", type=notExistingFile, required=True)
parser.add_argument("-s", "--resampleFactor",
                    help="Oversample the input image by this factor before processing",
                    type=int, default=5, required=False)
parser.add_argument("-t", "--threshold",
                    help="Values in the input image below this threshold will be considered as zero." +
                         " NB: the threshold is relative to the median of all the pixels in the image, " +
                         "so 0.1 means that all pixels less than ten percent of the median value are considered zero",
                    type=float, default=0.1, required=False)
parser.add_argument("-v", "--verbosity", required=False, default='info',
                    choices=['info', 'debug'])

parser.add_argument('--reverse', help="Use this to reverse the region selection (include instead of exclude)",
                    dest='reverse', action='store_true')
parser.add_argument('--no-reverse', help="This corresponds to the default choice (exclude regions)",
                    dest='reverse', action='store_false')
parser.set_defaults(reverse=True)

if __name__ == "__main__":
    args = parser.parse_args()

    log.setLevel(levels[args.verbosity])

    log.debug("Reading input image...")

    with pyfits.open(args.image) as f:

        data = f[0].data
        header = f[0].header

    log.debug("Input image shape: %s" % (list(data.shape)))

    medianExposure = numpy.median(data.flatten())
    threshold = args.threshold * medianExposure

    log.debug("Median value in input: %s" % (medianExposure))
    log.debug("threshold: %s" % (threshold))

    if (args.resampleFactor > 1):

        log.info("Resampling image by a factor of %s" % (args.resampleFactor))
        newData, newHeader = resampleImage(data, header,
                                           args.resampleFactor, threshold)
        log.info("New image shape: %s" % (list(newData.shape)))

    else:

        # Do not need resampling
        newData, newHeader = (data, header)

    pass

    log.debug("Reading region file...")
    regions = Regions(args.regionfile, newHeader, args.reverse)
    log.info("Found %s regions" % (regions.getNregions()))

    if regions.getNregions() > 0:

        # Filter
        log.debug("Getting mask...")
        mask = regions.getMask(newData.shape)

        log.debug("Filtering out pixels in regions...")
        filteredData = newData * (~mask)

    else:

        # No need to to any filtering
        filteredData = newData

    pass

    if (args.threshold != 0):
        log.debug("Zeroing values below threshold...")
        idx = filteredData <= threshold
        filteredData[idx] = 0.0

        if (args.makemask):
            log.debug("Masking...")
            filteredData[~idx] = 1

    log.debug("Writing output file...")
    pyfits.writeto(args.outfile, filteredData, newHeader)
