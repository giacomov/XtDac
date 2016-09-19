import numpy as np


class TimeInterval(object):
    def __init__(self, tstart, tstop, interestingRegion, rate_improvement, swap=False):

        self.tstart = tstart
        self.tstop = tstop
        self.interestingRegion = interestingRegion

        if (self.tstop <= self.tstart):

            if (swap):

                self.tstart = tstop
                self.tstop = tstart

            else:

                import pdb;
                pdb.set_trace()

                raise RuntimeError("Invalid time interval! TSTART must be before TSTOP and TSTOP-TSTART >0.")
        pass

        self.toBeRemoved = False

        # Note that the interesting region has still *all* the events within the region,
        # not just the one in this time interval, so we need to further filter down in time
        # to count the number of events in this time interval

        t = self.interestingRegion.box.getInRegionEventsTime()

        idx = (t >= self.tstart) & \
              (t <= self.tstop)

        self.nEvents = np.sum(idx)

        self.rateImprovement = float(rate_improvement)

        # box = interestingRegion.box

        # print("Evaluating region %s,%s for time interval %s-%s with probability %s" %(box.c1,box.c2,tstart,tstop,self.probability))

    pass

    def __repr__(self):

        return "Time interval: %s - %s (%s s)" % (self.tstart, self.tstop, self.tstop - self.tstart)

    def computeProbability(self):

        self.probability = self.interestingRegion.getProbability(self.tstart, self.tstop)

    def scheduleForRemoval(self):
        self.toBeRemoved = True

    def __eq__(self, other):
        if self.tstart == other.tstart and self.tstop == other.tstop and self.interestingRegion == other.interestingRegion:
            return True
        else:
            return False
        pass

    pass

    def getDuration(self):

        return self.tstop - self.tstart

    def getRate(self):

        return self.nEvents / float(self.getDuration())

    def getRateImprovement(self):

        return self.rateImprovement

    def reduceToIntersection(self, interval):

        if not self.overlapsWith(interval):

            self.tstart = None
            self.tstop = None

        else:

            self.tstart = max(self.tstart, interval.tstart)
            self.tstop = min(self.tstop, interval.tstop)

        pass

    pass

    def merge(self, interval):

        if (self.overlapsWith(interval)):

            self.tstart = min(self.tstart, interval.tstart)
            self.tstop = max(self.tstop, interval.tstop)

        else:

            raise RuntimeError("Could not merge non-overlapping intervals!")

    pass

    def overlapsWith(self, interval):

        if (interval.tstart == self.tstart or interval.tstop == self.tstop):
            return True

        if (interval.tstart > self.tstart and interval.tstart < self.tstop):
            return True

        if (interval.tstop > self.tstart and interval.tstop < self.tstop):
            return True

        if (interval.tstart < self.tstart and interval.tstop > self.tstop):
            return True

        return False


pass


class TimeIntervalConsolidator(object):
    def __init__(self, interestingRegions):

        self.interestingRegions = interestingRegions

    pass

    def consolidate(self, minNumberOfEvents=3):
        '''
    Keep only intervals above the floorProbability.
    Also, if there are overlapping time intervals from close-by regions,
    keep only the most significant one (the one with the lowest null-hyp.
    probability)
    '''

        # build a global list of intervals which survived the filter

        intervalsAboveFloorProb = []

        for reg in self.interestingRegions:

            thisEdges = reg.getIntervalsEdges()

            # Compute the rates (in the transformed frame) for this interesting region
            rates, tstarts, tstops = reg.getTransformedRates()

            # If there is only one interval, this cannot be an interesting region
            assert len(rates) >= 2

            # Compute rate improvements, i.e., the ratio between the rate in one interval and the rate in the
            # previous interval (or, if does not exist, in the following interval)

            rate_improvements = np.zeros_like(rates)
            rate_improvements[0] = rates[0] / rates[1]
            rate_improvements[1:] = rates[1:] / rates[:-1]

            # Make the list of TimeIntervals instances

            thisIntervals = []

            for rate_improvement, t1, t2, tt1, tt2 in zip(rate_improvements, thisEdges[:-1], thisEdges[1:],
                                                          tstarts, tstops):

                assert t1==tt1 and t2==tt2

                thisIntervals.append(TimeInterval(t1, t2, reg, rate_improvement))

            # Filter by probability
            # thisIntervalsClean = filter(lambda interval: interval.probability < floorProbability, thisIntervals)

            # Filter by minimum number of events
            thisIntervalsClean = filter(lambda interval:
                                        interval.nEvents > minNumberOfEvents and interval.getRateImprovement() >= 1.0,
                                        thisIntervals)

            intervalsAboveFloorProb.extend(thisIntervalsClean)
        pass

        # Now check if there are overlapping intervals
        # (remember: there is a flag in the TimeInterval class
        # called "toBeRemoved" which we will change if needed
        # in the following loop)

        for int1 in intervalsAboveFloorProb:

            # The following can happen if we are not in the first
            # iteration, and the current interval has already been
            # flagged for removal. Don't need to process it again
            if (int1.toBeRemoved):
                continue

            for int2 in intervalsAboveFloorProb:

                if (int1 == int2):
                    continue

                # Check if interval1 and interval2 overlaps in time

                if int1.overlapsWith(int2):

                    # Check if int1 and int2 are from close-by regions

                    if (int1.interestingRegion.overlapsWith(int2.interestingRegion)):

                        # Compute their probabilities
                        # int1.computeProbability()
                        # int2.computeProbability()

                        # Yes, they are from overlapping regions. Choose the one with the lowest probability
                        if int1.getRateImprovement() > int2.getRateImprovement():

                            int2.toBeRemoved = True

                        else:

                            # Remove int1
                            int1.toBeRemoved = True

        nonOverlappingIntervals = filter(lambda x: x.toBeRemoved == False, intervalsAboveFloorProb)

        return nonOverlappingIntervals

    pass
