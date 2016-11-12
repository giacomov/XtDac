import math
import sys

import numpy
import scipy.stats

import Likelihood
import cellDetect

try:

    import astropy.io.fits as pyfits

except:

    import pyfits

import logging

log = logging.getLogger("main")


class FixedBinSearch(object):

    def __init__(self, X, Y, t, timebin, eventfile, expomap, hwu, region=None):

        # Store the hardware unit
        self.hwu = hwu

        # Read the events
        self.X = X
        self.Y = Y
        self.time = t

        self.timebin = float(timebin)

        # Read the source catalog
        sources = []
        with pyfits.open(eventfile) as f:

            srcX = []
            srcY = []

            for row in f['REGION'].data:
                srcX.append(row.field("X")[0])
                srcY.append(row.field("Y")[0])

        sources = numpy.array(zip(srcX, srcY))

        # Transform RA,Dec in X,Y for the sources

        self.eventfile = eventfile

        # wcs = XMMWCS.XMMWCS(self.eventfile)

        # sources = numpy.array( wcs.sky2xy(zip(ras,decs)) )

        log.debug("Read %s sources from the source catalog" % (len(sources)))

        if region is None:

            self.xmin = self.X.min()
            self.xmax = self.X.max()

            self.ymin = self.Y.min()
            self.ymax = self.Y.max()

            # Now keep only the sources within this region
            idx = ((sources[:, 0] >= self.xmin) &
                   (sources[:, 0] <= self.xmax) &
                   (sources[:, 1] >= self.ymin) &
                   (sources[:, 1] <= self.ymax))

            self.sources = sources[idx]

        else:

            # Use the region to filter the sources

            filt = region.get_filter()

            self.sources = []

            for src in sources:

                if filt.inside1(src[0], src[1]):
                    self.sources.append(src)

            # Use the region to determine the extreme of X and Y

            shape = region[0]

            # NOTE: Assuming all units are in pixels

            if shape.name == 'circle':

                xc, yc, radius = shape.coord_list

                # sometimes for no reason the radius is negative. Fix that
                radius = abs(radius)

            elif shape.name == 'ellipse':

                xc, yc, r1, r2, ang = shape.coord_list

                radius = max(abs(r1), abs(r2))

            else:

                raise RuntimeError("Region of shape %s is not supported" % (shape.name))

            pass

            self.xmin = xc - radius
            self.xmax = xc + radius

            self.ymin = yc - radius
            self.ymax = yc + radius

        nsrc = len(self.sources)

        nmax = 100

        if nsrc > nmax:
            log.info("There were %s sources in the region, they are too many. Keeping only %s" % (nsrc, nmax))

            self.sources = self.sources[:nmax]

        log.debug("Kept %s source(s) within the region" % (len(self.sources)))

        self.expomap = expomap

    def scan(self, regfilt, threshold=5.0):

        bins = numpy.arange(self.time.min(), self.time.max() + 1e-4, self.timebin)

        xbins = numpy.arange(self.xmin, self.xmax, 25)
        ybins = numpy.arange(self.ymin, self.ymax, 25)

        region = Likelihood.Region(self.xmin, self.xmax,
                                   self.ymin, self.ymax,
                                   40, regfilt,
                                   self.expomap, self.eventfile)

        ls = Likelihood.Likelihood(self.X, self.Y, region)

        # Build the likelihood model

        # Bkg model (isotropic)
        iso = Likelihood.Isotropic("bkg", 1.0)

        m = Likelihood.GlobalModel("likeModel") + iso

        knownSources = []

        for i, src in enumerate(self.sources):

            thisSrc = Likelihood.PointSource("src%i" % (i + 1), src[0], src[1], self.eventfile, self.hwu)

            m += thisSrc

            if i >= 10:
                log.info("Fixing normalization of source %s" % (i + 1))

                thisSrc["src%i_norm" % (i + 1)].fix()

            knownSources.append(thisSrc)

        ls.setModel(m)

        # Minimize the mlogLike to get the background level
        like0 = ls.minimize()

        figIntegrated = ls.plot()

        excessesFound = []

        figs = [figIntegrated]

        for t1, t2 in zip(bins[:-1], bins[1:]):

            sys.stderr.write(".")

            idx = (self.time >= t1) & (self.time < t2)

            # db = DbScanner(self.X[idx],self.Y[idx])

            # db = SAODetect(self.X[idx],self.Y[idx])

            # db = cellDetect.CellDetect(self.X[idx],self.Y[idx],200,3)

            ls = Likelihood.Likelihood(self.X[idx], self.Y[idx], region)

            totCounts = ls.getTotalCounts()

            if totCounts < 3:
                # No need to perform any fit. Too few counts in this
                # interval

                log.debug("Skip interval %s - %s, less than 3 counts here" % (t1, t2))

                continue

            # Build the likelihood model

            # Bkg model (isotropic)
            iso = Likelihood.Isotropic("bkg", 0.1)

            m = Likelihood.GlobalModel("likeModel") + iso

            for i, src in enumerate(knownSources):
                m += src

            ls.setModel(m)

            # Minimize the mlogLike to get the background level
            like0 = ls.minimize(verbose=0)

            # ls.plot()

            bkgLevel = iso['bkg_amp'].value / pow(region.binsize, 2.0)
            error = iso['bkg_amp'].error / pow(region.binsize, 2.0)

            # print("Bkg: %s +/- %s" %(bkgLevel, error))

            db = cellDetect.CellDetect(self.X[idx], self.xmin, self.xmax,
                                       self.Y[idx], self.ymin, self.ymax,
                                       200, 3)

            centroids = db.findExcesses(bkgLevel * pow(200, 2), error * pow(200, 2))

            newSources = []

            if len(centroids) > 0:
                # Verify if this source is truly significant

                TSs = []

                for (x, y) in centroids:

                    sys.stderr.write("-")

                    # Avoid computing the TS for a source too close to the ones already
                    # known

                    d = map(lambda src: math.sqrt(pow(src.xpos - x, 2) + pow(src.ypos - y, 2)), knownSources)

                    dd = filter(lambda x: x < 5.0 / 0.05, d)

                    if len(dd) > 0:
                        # There is a source close to this centroid. No need
                        # to investigate further
                        continue

                    thisSrc = Likelihood.PointSource("testsrc", x, y, self.eventfile, self.hwu,
                                                     knownSources[0].outfile)

                    ls.model += thisSrc

                    like1 = ls.minimize(verbose=0)

                    TSs.append(2 * (like0 - like1))
                    # print("(X,Y) = (%s,%s) -> TS = %s" %(x,y,TSs[-1]))

                    ls.model.removeSource("testsrc")

                    # Using the approximation TS = sqrt( significance )

                    if TSs[-1] >= pow(threshold, 2):
                        newSources.append([x, y, TSs[-1]])

            if len(newSources) > 0:
                db.centroids = numpy.array(newSources)
                fig = db.plot(xbins, ybins)

                figs.append(fig)

                for (x, y, ts) in newSources:
                    excessesFound.append([t1, t2, x, y, math.sqrt(ts)])

        sys.stderr.write("\n")

        return excessesFound, figs


def get_pvalue(significance):
    return scipy.stats.norm().sf(significance)


def consolidateTimeIntervals(allExcesses, pixelSize):
    h = 0

    while True:

        if h >= allExcesses.shape[0]:
            break

        # Take h excess as reference
        ref = allExcesses[h]

        # Compute the distance between the reference
        # and all the others
        d = numpy.sqrt(numpy.power(ref[2] - allExcesses[h:, 2], 2) +
                       numpy.power(ref[3] - allExcesses[h:, 3], 2))

        idx = (d <= 10.0 / pixelSize)

        if numpy.sum(idx) > 0:

            # Select all excesses close to the reference

            thisExcesses = allExcesses[h:][idx]

            # If they are contiguous also in the time domain,
            # collapse them

            toBeRemoved = []

            for i, exc in enumerate(thisExcesses):

                if exc[0] - ref[1] == 0:
                    # Contiguous also in time
                    # Collapse these excesses to one

                    # Update end time

                    ref[1] = exc[1]

                    # Update significance by combining the two significances
                    # (formula from:
                    # http://www.loujost.com/Statistics%20and%20Physics/Significance%20Levels/CombiningPValues.htm)

                    # If the sigma value in exc[-1] is too high, the get_pvalue function
                    # will return zero because the probability is so low that it underflows
                    # the precision of the macine. Hence, we put a lower limit at 1e-20

                    k1 = max(get_pvalue(exc[-1]), 1e-20)

                    k0 = max(get_pvalue(ref[-1]), 1e-20)

                    k = k0 * k1

                    new_pvalue = k - k * math.log(k)

                    ref[-1] = scipy.stats.norm().isf(new_pvalue)

                    toBeRemoved.append(i + h)

            allExcesses = numpy.delete(allExcesses, toBeRemoved, axis=0)

        h += 1

    return allExcesses
