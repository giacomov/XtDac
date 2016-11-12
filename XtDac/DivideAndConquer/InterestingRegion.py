import numpy
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import XMMWCS


class InterestingRegion(object):
    def __init__(self, box, intervals, time, X, Y, tstart, tstop):

        self.box = box
        self.box.readEvents(True, time, X, Y, tstart, tstop)

        self.intervals = intervals

    pass

    def overlapsWith(self, intreg2):

        # Assume there is no rotation
        assert (self.box.rot == 0)

        X1, Y1, W1, H1 = (self.box.c1, self.box.c2, self.box.width, self.box.height)
        X2, Y2, W2, H2 = (intreg2.box.c1, intreg2.box.c2, intreg2.box.width, intreg2.box.height)

        if (X1 <= X2 + W2 and
                        X1 + W1 >= X2 and
                    Y1 <= Y2 + H2 and
                        Y1 + H1 >= Y2):

            # Overlap
            return True

        else:

            return False

        pass

    pass

    def getIntervalsEdges(self):
        return self.intervals

    def getWidthAndHeight(self):
        return (self.box.width, self.box.height)

    def getCenterSkyCoordinates(self):

        wcs = XMMWCS.XMMWCS(self.box.eventFile)
        ra, dec = wcs.xy2sky([self.getCenterPhysicalCoordinates()])[0]

        return ra, dec

    pass

    def getCenterPhysicalCoordinates(self):

        return [self.box.c1 + self.box.width / 2.0,
                self.box.c2 + self.box.height / 2.0]

    def _computeLightCurve(self):

        if (hasattr(self, "lightCurveComputed")):
            # Don't need to do it again
            return
        pass

        # obtain integral distribution
        self.intDistr = self.box.buildBkgIntegralDistribution()

        # Transform the arrival times of the events to a time system where
        # the background is flat
        tt = self.intDistr(self.box.getInRegionEventsTime())
        intervals_ = self.intDistr(self.intervals)
        widths_ = intervals_[1:] - intervals_[:-1]

        # Make the light curve in the transformed frame
        (counts, _) = numpy.histogram(tt, intervals_)

        self.rates = counts / widths_

        # Plot the light curve with the time in the original frame
        self.leftEdges = self.intervals[:-1]
        self.rightEdges = self.intervals[1:]
        self.widths = self.rightEdges - self.leftEdges

        self.lightCurveComputed = True

    def getTransformedRates(self):
        """
        Return the rate of events in the transformed frame. Note that the absolute value doesn't mean anything,
        only the relative value is meaningful to measure an increase with respect to the background

        :return: tstart, tstop, rates (tstart and tstop are in the original frame)
        """

        self._computeLightCurve()

        return self.rates, self.leftEdges, self.rightEdges

    def getLightCurve(self, **kwargs):
        self._computeLightCurve()

        figure = plt.figure(**kwargs)
        sub = figure.add_subplot(111)
        sub.bar(self.leftEdges, self.rates, self.widths, color='white', edgecolor='black')
        sub.set_xlabel("Mission Elapsed Time (s)")
        sub.set_ylabel("Rate")
        sub.set_title("Detrended light curve")
        return figure

    pass

    # def getProbability(self, t1, t2):
    #     self._computeLightCurve()
    #
    #     p = self.box.getProbability(t1, t2, integraldistribution=self.intDistr)
    #
    #     return p
    #
    # pass

    def getStamps(self, **kwargs):

        # For each interval, create image

        inRegionEventsTime = self.box.getInRegionEventsTime()
        inRegionEventsX = self.box.getInRegionEventsX()
        inRegionEventsY = self.box.getInRegionEventsY()

        # Decide binning extremes and re-binning factor
        (xmin, xmax) = (self.box.c1, self.box.c1 + self.box.width)
        (ymin, ymax) = (self.box.c2, self.box.c2 + self.box.height)

        # Try to get a reasonable number of bins based on the size of the region, but
        # do never do less than 20 bins

        nbinsX = max(numpy.rint((xmax - xmin) / 40), 20)
        nbinsY = max(numpy.rint((ymax - ymin) / 40), 20)

        # Decide how many rows and columns to have in the figure
        l = len(self.intervals)
        nrows = 0
        ncols = 0
        for i in range(l):
            if (nrows * ncols >= l):
                break
            nrows += 1
            if (nrows * ncols >= l):
                break
            ncols += 1

        fig, subs = plt.subplots(nrows, ncols, **kwargs)

        for sub, t1, t2 in zip(subs.flatten(), self.intervals[:-1], self.intervals[1:]):
            thisEventsIdx = (inRegionEventsTime >= t1) & (inRegionEventsTime < t2)
            thisTimes = inRegionEventsTime[thisEventsIdx]
            thisX = inRegionEventsX[thisEventsIdx]
            thisY = inRegionEventsY[thisEventsIdx]

            sub.hist2d(thisX, thisY,
                       [nbinsX, nbinsY],
                       range=[[xmin, xmax], [ymin, ymax]],
                       cmap='afmhot')
        pass

        # Remove all axis from the subplots
        for sub in subs.flatten():
            sub.xaxis.set_visible(False)
            sub.yaxis.set_visible(False)
        pass

        return fig

    def getFindingMap(self, **kwargs):
        # Plot a figure with the whole event file and the position of the region

        fig = plt.figure(**kwargs)
        sub = fig.add_subplot(111)
        sub.hist2d(self.box.x, self.box.y, [200, 200], cmap='afmhot')
        # Now overlopt the region
        rect = Rectangle((self.box.c1, self.box.c2),
                         self.box.width, self.box.height,
                         color='green', linewidth=2.0,
                         alpha=0.5)
        sub.add_patch(rect)
        sub.xaxis.set_visible(False)
        sub.yaxis.set_visible(False)

        return fig

    def clearMemory(self):
        self.box.clearMemory()


pass
