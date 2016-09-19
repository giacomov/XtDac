# Author: Giacomo Vianello (giacomov@stanford.edu)

import os
import sys

import numpy
import numexpr

try:
    import scipy.stats
except:
    pass

try:
    import astropy.io.fits as pyfits
except:
    # If this fail there is no way out
    import pyfits
pass

from xtwp4.BayesianBlocks import bayesian_blocks, bayesian_blocks_not_unique


class Box(object):

    def __init__(self, c1, c2, width, height, rot, coordsys):
        '''
        A Box region. It can be in fk5, in which case
        height and width are in arcsec and c1 and c2 are respectively R.A.,Dec;
        or in physical coordinates, in which case height and width are in pixels,
        and c1 and c2 in pixels as well (detector coordinates). c1 and c2 are the
        coordinates of the lower left corner of the box.
        '''
        self.c1 = c1
        self.c2 = c2
        self.width = width
        self.height = height
        self.rot = rot

        # This will become true when events are selected and kept as members
        # of the class. Is set back on False by the clearMemory method
        self.eventsLoaded = False

        self.id = ''

        self.coordsys = coordsys
        allowed = ['physical', 'sky']
        if (coordsys not in allowed):
            raise ValueError("coordsys '%s' is unrecognized. Legal values are: %s." % (coordsys, ",".join(allowed)))

    pass

    def setId(self, regid):
        self.id = regid

    pass

    def setEventfile(self, eventFile):

        if (not os.path.exists(eventFile)):
            raise IOError("File %s does not exist or is not accesible." % (eventFile))

        else:

            self.eventFile = eventFile

    pass

    def readEvents(self, preLoaded=False, time=None, X=None, Y=None, tstart=None, tstop=None):

        if (self.coordsys == "physical"):
            xcoord = "DETX"
            ycoord = "DETY"
        elif (self.coordsys == "sky"):
            xcoord = "X"
            ycoord = "Y"
        pass

        if (preLoaded == False):
            # Read the file
            with pyfits.open(self.eventFile) as fitsFile:
                sys.stderr.write("\nReading FITS file...\n")
                # Select the events
                X = fitsFile['EVENTS'].data.field(xcoord)
                Y = fitsFile['EVENTS'].data.field(ycoord)
                time = fitsFile['EVENTS'].data.field('TIME')

                # Read start and stop time of the observation
                self.tstart = fitsFile['EVENTS'].header.get("TSTART")
                self.tstop = fitsFile['EVENTS'].header.get("TSTOP")

                self.eventsLoaded = True

        else:
            self.tstart = tstart
            self.tstop = tstop
            self.eventsLoaded = True
        pass

        # Using numexpr.eval because it is much faster than the normal
        # pyfits/numpy way. Indeed, normally an expression with many ANDs
        # is evaluated in pieces and then ANDed, while numexpr evaluate the
        # whole expression once for each entry
        expr = '( (X >= %s) & (X < %s) & (Y >= %s) & (Y < %s))' % (
            self.c1, self.c1 + self.width,
            self.c2, self.c2 + self.height)

        idx = numexpr.evaluate(expr)

        self.arrivalTimes = time
        self.x = X
        self.y = Y
        self.inRegionIdx = idx
        self.nInRegion = numpy.sum(self.inRegionIdx)

        # Use all the events outside the region as background
        # (although only a subset might be used in the end, see below)
        # NB: ~ means "logical not"

        self.setBackground(~idx)

        # self.outOfRegionIdx     = ~idx #Logical not
        # self.nOutOfRegion       = numpy.sum(self.outOfRegionIdx)

        if (self.nInRegion == 0):

            self.empty = True
            self.filledArea = 0

        else:

            self.empty = False

            # If we have many background events, choose among them
            # the one closest to the region of interest

            if (self.nOutOfRegion > 20 * self.nInRegion):
                # We have more than 20 times the events in the region which can be used for background
                # estimation. Therefore, we have the luxury of choosing the events.
                # Let's use the events closest to this region, which better reflect the local background

                # Compute the center of this region
                xc, yc = (self.c1 + self.width / 2.0, self.c2 + self.height / 2.0)

                # Do not take the sqrt for speed
                c1 = self.c1
                c2 = self.c2
                x = self.x
                y = self.y
                distances = numexpr.evaluate('''(xc - x)**2 + (yc - y)**2''')

                # Of all the events outside the region, select the ones closest to the center
                idx = numpy.argsort(distances)
                howMany = 20 * self.nInRegion
                not_chosen = idx[howMany:]
                self.outOfRegionIdx[not_chosen] = False
                self.nOutOfRegion = numpy.sum(self.outOfRegionIdx)
            pass

            if (self.x.shape[0] < 50):

                # we don't have enought events to estimate the filled area
                self.filledArea = 1

            else:

                # Rought computation of filled area: regions on the border of the
                # image will have a low filled area
                evts_minX = self.x[self.inRegionIdx].min()
                evts_maxX = self.x[self.inRegionIdx].max()
                evts_minY = self.y[self.inRegionIdx].min()
                evts_maxY = self.y[self.inRegionIdx].max()

                area = (evts_maxX - evts_minX) * (evts_maxY - evts_minY)

                self.filledArea = area / (self.width * self.height)
        pass

    pass

    def setBackground(self, idx):
        '''
        This overwrite the background selection made by the readEvents method
        '''

        self.outOfRegionIdx = idx
        self.nOutOfRegion = numpy.sum(self.outOfRegionIdx)

    pass

    def clearMemory(self):

        del self.arrivalTimes
        del self.x
        del self.y
        del self.inRegionIdx
        del self.outOfRegionIdx

        self.eventsLoaded = False

    def getInRegionEventsTime(self):
        if (not self.eventsLoaded):
            self.readEvents()
        pass
        return self.arrivalTimes[self.inRegionIdx]

    def getInRegionEventsX(self):
        if (not self.eventsLoaded):
            self.readEvents()
        pass
        return self.x[self.inRegionIdx]

    def getInRegionEventsY(self):
        if (not self.eventsLoaded):
            self.readEvents()
        pass
        return self.y[self.inRegionIdx]

    def isEmpty(self):
        return self.empty

    def buildBkgIntegralDistribution(self, unbinned=True):
        # Make the integral distribution of the background
        if (not self.eventsLoaded):
            self.readEvents()
        pass

        if (unbinned):

            # Add a point at the beginning and at the end to be sure that the
            # integral distribution will be extrapolated correctly if needed
            x = numpy.copy(self.arrivalTimes[self.outOfRegionIdx])

            x = numpy.insert(x,
                             [0, x.shape[0]],
                             [self.tstart,
                              self.tstop])

            x.sort()

            cumulative = numpy.arange(1.0, self.nOutOfRegion + 1, 1)
            y = numpy.insert(cumulative, [0, cumulative.shape[0]], [0, self.nOutOfRegion + 1])
            y.sort()

        else:

            # Make a light curve with 25 events in each bin
            edges = self.arrivalTimes[self.outOfRegionIdx][::25][1:]
            edges = numpy.insert(edges, [0, edges.shape[0]], [self.tstart, self.tstop])

            res = numpy.histogram(self.arrivalTimes[self.outOfRegionIdx], edges)

            x = ((edges[1:] + edges[:-1]) / 2.0)
            x = numpy.insert(x, [0, x.shape[0]], [self.tstart, self.tstop])

            y = numpy.cumsum(res[0])
            y = numpy.insert(y, [0, y.shape[0]], [0, self.nOutOfRegion + 1])

        # Build the interpolating function for the integral distribution
        # intDist(x) will give the number of counts between 0 and x

        # Cannot save intDistr as attribute of the class because it would not allow
        # to pickle the class, which will make it unusable with the multiprocessing
        # python facilities

        # intDistr              = scipy.interpolate.interp1d(x,y)
        intDistr = lambda xx: numpy.interp(xx, x, y)

        return intDistr

    def findExcesses(self, nullHypProb):

        '''Run the Bayesian Block algorithm and returns the list of intervals'''

        intDistr = self.buildBkgIntegralDistribution()
        tt = self.arrivalTimes[self.inRegionIdx]

        if tt.shape[0] != numpy.unique(tt).shape[0]:

            raise RuntimeError("Non-unique time stamps in event file")

            # res = bayesian_blocks_not_unique(tt, self.tstart, self.tstop, nullHypProb)

        else:

            res = bayesian_blocks(tt, self.tstart, self.tstop, nullHypProb, intDistr)

        return res

    pass

    def poissonProbability(self, nobs, npred):

        # The probability of obtaining nobs *or more* when npred is expected

        return (scipy.stats.distributions.poisson.sf(nobs, npred) +
                scipy.stats.distributions.poisson.pmf(nobs, npred))

    def getProbability(self, tstart, tstop, **kwargs):

        intDistr = None
        onlyExcesses = True

        for k, v in kwargs.iteritems():

            if (k.lower() == "integraldistribution"):

                intDistr = v

            elif (k.lower() == "onlyexcesses"):

                onlyExcesses = bool(v)

            pass

        pass

        if (intDistr == None):
            intDistr = self.buildBkgIntegralDistribution()

        # How many out-of-region events outside this interval
        tt = self.arrivalTimes[self.outOfRegionIdx]
        idx = numexpr.evaluate("(tt < tstart) | (tt > tstop)")
        nout = numpy.sum(idx)

        # How many in-region events outside this interval
        tt = self.arrivalTimes[self.inRegionIdx]
        idx = numexpr.evaluate("(tt < tstart) | (tt > tstop)")
        nin = numpy.sum(idx)

        approxRatio = float(nout) / float(nin)

        # if approxRatio is zero and there are more than 3 counts as nin,
        # use 1e-7 as probability (arbitrary value) just to keep this interval
        # alive for further screening

        if approxRatio == 0:

            if nin >= 3:

                p = 1.01e-7

            else:

                p = 1

        else:

            Npred = (intDistr(tstop) - intDistr(tstart)) / approxRatio

            tt = self.arrivalTimes[self.inRegionIdx]
            idx = numexpr.evaluate("(tt >= tstart) & (tt <= tstop)")
            Nobs = numpy.sum(idx)

            if (onlyExcesses and Nobs <= Npred):

                p = 1

            else:

                p = self.poissonProbability(Nobs, Npred)

        return p

    pass

    def __repr__(self):
        return 'Upper corner: (%s,%s), width: %s, height: %s, coord.sys: %s' % (
        self.c1, self.c2, self.width, self.height, self.coordsys)

    pass

    def _getXYZ(self, th, ph):
        theta = numpy.deg2rad(th)
        phi = numpy.deg2rad(ph)
        x = numpy.sin(phi) * numpy.cos(theta)
        y = numpy.sin(phi) * numpy.sin(theta)
        z = numpy.cos(phi)
        return x, y, z

    pass

    def _getThetaPhi(self, x, y, z):
        theta = numpy.arccos(z)
        phi = numpy.arctan2(y, x)
        return numpy.rad2deg(phi), numpy.rad2deg(theta)

    pass

    def rotate(self, xc, yc, angle):
        '''
        Rotate the system around the point (xc,yc) of a certain angle (in deg).
        '''

        # Translate to the new origin
        newC1 = self.c1 - xc
        newC2 = self.c2 - yc

        # Rotate
        angleRad = numpy.deg2rad(angle)
        newC1 = newC1 * numpy.cos(angleRad) - newC2 * numpy.sin(angleRad)
        newC2 = newC1 * numpy.sin(angleRad) + newC2 * numpy.cos(angleRad)

        # Translate back to the old origin
        self.c1 = newC1 + xc
        self.c2 = newC2 + yc
        self.rot = angle

    pass

    def rotateSky(self, ra_c, dec_c, deg):
        '''
        Rotate the box of 'deg' degrees around the point with R.A., Dec = (ra_c,dec_c).
        Only works for fk5 coordinates.
        '''
        if (self.coordsys == 'physical'):
            raise RuntimeError("Cannot rotate physical coordinates")
        pass

        x_c, y_c, z_c = self._getXYZ(ra_c, dec_c)
        x, y, z = self._getXYZ(self.c1, self.c2)

        rads = numpy.deg2rad(deg)

        cosa = numpy.cos(rads)
        sina = numpy.sin(rads)

        x1 = (x * cosa + (1 - cosa) * (x_c * x_c * x + x_c * y_c * y + x_c * z_c * z)
              + (y_c * z - z_c * y) * sina)
        y1 = (y * cosa + (1 - cosa) * (y_c * x_c * x + y_c * y_c * y + y_c * z_c * z)
              + (z_c * x - x_c * z) * sina)
        z1 = (z * cosa + (1 - cosa) * (z_c * x_c * x + z_c * y_c * y + z_c * z_c * z)
              + (x_c * y - y_c * x) * sina)

        ra1, dec1 = self._getThetaPhi(x1, y1, z1)
        self.c1 = ra1
        self.c2 = dec1

    pass

    def getDs9Region(self):
        # Set the units for the width and height for the region file
        units = ''
        conv = 1

        string = []
        # Ds9 boxes are defined starting from the center of the box
        string.append("physical\n")
        string.append('''box(%s,%s,%s%s,%s%s,%s)\n''' % ((self.c1 + self.width / 2.0 / conv),
                                                         (self.c2 + self.height / 2.0 / conv),
                                                         self.width, units, self.height,
                                                         units, self.rot))
        return string

    pass

    def writeRegion(self, outfile):
        with open(outfile, 'w+') as f:
            f.write("\n".join(self.getDs9Region()))

    pass

    def writeFitsFile(self, tstart, tstop, outfile):
        idx = (self.arrivalTimes >= tstart) & (self.arrivalTimes < tstop) & self.inRegionIdx
        xs = self.x[idx]
        ys = self.y[idx]
        ts = self.arrivalTimes[idx]

        col1 = pyfits.Column(name='X', format='E', array=xs)
        col2 = pyfits.Column(name='Y', format='E', array=ys)
        col3 = pyfits.Column(name='TIME', format='E', array=ts)
        coldefs = pyfits.ColDefs([col1, col2, col3])
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
        tbhdu.header.set("EXTNAME", "EVENTS")

        tbhdu.writeto(outfile, clobber=True)


pass
