# Author: Giacomo Vianello (giacomov@stanford.edu)
# Copyright 2014 EXTraS.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#   http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import sys

import numpy as np
import numexpr

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("bayesian_blocks")

__all__ = ['bayesian_blocks', 'bayesian_blocks_not_unique']


def bayesian_blocks_not_unique(tt, ttstart, ttstop, p0):

    # Verify that the input array is one-dimensional
    tt = np.asarray(tt, dtype=float)

    assert tt.ndim == 1

    # Now create the array of unique times

    unique_t = np.unique(tt)

    t = tt
    tstart = ttstart
    tstop = ttstop

    # Create initial cell edges (Voronoi tessellation) using the unique time stamps

    edges = np.concatenate([[tstart],
                            0.5 * (unique_t[1:] + unique_t[:-1]),
                            [tstop]])

    # The last block length is 0 by definition
    block_length = tstop - edges

    if np.sum((block_length <= 0)) > 1:

        raise RuntimeError("Events appears to be out of order! Check for order, or duplicated events.")

    N = unique_t.shape[0]

    # arrays to store the best configuration
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)

    # Pre-computed priors (for speed)
    # eq. 21 from Scargle 2012

    priors = 4 - np.log(73.53 * p0 * np.power(np.arange(1, N + 1), -0.478))

    # Count how many events are in each Voronoi cell

    x, _ = np.histogram(t, edges)

    # Speed tricks: resolve once for all the functions which will be used
    # in the loop
    cumsum = np.cumsum
    log = np.log
    argmax = np.argmax
    numexpr_evaluate = numexpr.evaluate
    arange = np.arange

    # Decide the step for reporting progress
    incr = max(int(float(N) / 100.0 * 10), 1)

    logger.debug("Finding blocks...")

    # This is where the computation happens. Following Scargle et al. 2012.
    # This loop has been optimized for speed:
    # * the expression for the fitness function has been rewritten to
    #  avoid multiple log computations, and to avoid power computations
    # * the use of scipy.weave and numexpr has been evaluated. The latter
    #  gives a big gain (~40%) if used for the fitness function. No other
    #  gain is obtained by using it anywhere else

    # Set numexpr precision to low (more than enough for us), which is
    # faster than high
    oldaccuracy = numexpr.set_vml_accuracy_mode('low')
    numexpr.set_num_threads(1)
    numexpr.set_vml_num_threads(1)

    for R in range(N):
        br = block_length[R + 1]
        T_k = block_length[:R + 1] - br

        # N_k: number of elements in each block
        # This expression has been simplified for the case of
        # unbinned events (i.e., one element in each block)
        # It was:
        N_k = cumsum(x[:R + 1][::-1])[::-1]
        # Now it is:
        #N_k = arange(R + 1, 0, -1)

        # Evaluate fitness function
        # This is the slowest part, which I'm speeding up by using
        # numexpr. It provides a ~40% gain in execution speed.

        fit_vec = numexpr_evaluate('''N_k * log(N_k/ T_k) ''',
                                   optimization='aggressive',
                                   local_dict={'N_k': N_k, 'T_k': T_k})

        p = priors[R]

        A_R = fit_vec - p

        A_R[1:] += best[:R]

        i_max = argmax(A_R)

        last[R] = i_max
        best[R] = A_R[i_max]

    pass

    numexpr.set_vml_accuracy_mode(oldaccuracy)

    logger.debug("Done\n")

    # Now find blocks
    change_points = np.zeros(N, dtype=int)
    i_cp = N
    ind = N
    while True:
        i_cp -= 1
        change_points[i_cp] = ind

        if ind == 0:
            break

        ind = last[ind - 1]

    change_points = change_points[i_cp:]

    finalEdges = edges[change_points]

    return np.asarray(finalEdges)

def bayesian_blocks(tt, ttstart, ttstop, p0, bkgIntegralDistr=None, myLikelihood=None):
    """Divide a series of events characterized by their arrival time in blocks
    of perceptibly constant count rate. If the background integral distribution
    is given, divide the series in blocks where the difference with respect to
    the background is perceptibly constant.


    Args:
      tt (iterable): An iterable (list, numpy.array...) containing the arrival
                     time of the events.
                     NOTE: the input array MUST be time-ordered, and without
                     duplicated entries. To ensure this, you may execute the
                     following code:

                     tt_array = numpy.asarray(tt)
                     tt_array = numpy.unique(tt_array)
                     tt_array.sort()

                     before running the algorithm.

      p0 (float): The probability of finding a variations (i.e., creating a new
                  block) when there is none. In other words, the probability of
                  a Type I error, i.e., rejecting the null-hypothesis when is
                  true. All found variations will have a post-trial significance
                  larger than p0.

      bkgIntegralDistr (function, optional): the integral distribution for the
                  background counts. It must be a function of the form f(x),
                  which must return the integral number of counts expected from
                  the background component between time 0 and x.

    Returns:
      numpy.array: the edges of the blocks found

    """

    # Verify that the input array is one-dimensional
    tt = np.asarray(tt, dtype=float)

    assert tt.ndim == 1

    if (bkgIntegralDistr is not None):
        # Transforming the inhomogeneous Poisson process into an homogeneous one with rate 1,
        # by changing the time axis according to the background rate
        logger.debug("Transforming the inhomogeneous Poisson process to a homogeneous one with rate 1...")
        t = np.array(bkgIntegralDistr(tt))
        logger.debug("done")

        # Now compute the start and stop time in the new system
        tstart = bkgIntegralDistr(ttstart)
        tstop = bkgIntegralDistr(ttstop)
    else:
        t = tt
        tstart = ttstart
        tstop = ttstop
    pass

    # Create initial cell edges (Voronoi tessellation)
    edges = np.concatenate([[tstart],
                            0.5 * (t[1:] + t[:-1]),
                            [tstop]])

    # Create the edges also in the original time system
    edges_ = np.concatenate([[ttstart],
                             0.5 * (tt[1:] + tt[:-1]),
                             [ttstop]])


    # Create a lookup table to be able to transform back from the transformed system
    # to the original one
    lookupTable = {key: value for (key, value) in zip(edges, edges_)}

    # The last block length is 0 by definition
    block_length = tstop - edges

    if np.sum((block_length <= 0)) > 1:

        raise RuntimeError("Events appears to be out of order! Check for order, or duplicated events.")

    N = t.shape[0]

    # arrays to store the best configuration
    best = np.zeros(N, dtype=float)
    last = np.zeros(N, dtype=int)
    best_new = np.zeros(N, dtype=float)
    last_new = np.zeros(N, dtype=int)

    # Pre-computed priors (for speed)

    if (myLikelihood):

        priors = myLikelihood.getPriors(N, p0)

    else:

        # eq. 21 from Scargle 2012
        priors = 4 - np.log(73.53 * p0 * np.power(np.arange(1, N + 1), -0.478))
    pass

    x = np.ones(N)

    # Speed tricks: resolve once for all the functions which will be used
    # in the loop
    cumsum = np.cumsum
    log = np.log
    argmax = np.argmax
    numexpr_evaluate = numexpr.evaluate
    arange = np.arange

    # Decide the step for reporting progress
    incr = max(int(float(N) / 100.0 * 10), 1)

    logger.debug("Finding blocks...")

    # This is where the computation happens. Following Scargle et al. 2012.
    # This loop has been optimized for speed:
    # * the expression for the fitness function has been rewritten to
    #  avoid multiple log computations, and to avoid power computations
    # * the use of scipy.weave and numexpr has been evaluated. The latter
    #  gives a big gain (~40%) if used for the fitness function. No other
    #  gain is obtained by using it anywhere else

    times = []
    TSs = []

    # Set numexpr precision to low (more than enough for us), which is
    # faster than high
    oldaccuracy = numexpr.set_vml_accuracy_mode('low')
    numexpr.set_num_threads(1)
    numexpr.set_vml_num_threads(1)

    for R in range(N):
        br = block_length[R + 1]
        T_k = block_length[:R + 1] - br

        # N_k: number of elements in each block
        # This expression has been simplified for the case of
        # unbinned events (i.e., one element in each block)
        # It was:
        # N_k = cumsum(x[:R + 1][::-1])[::-1]
        # Now it is:
        N_k = arange(R + 1, 0, -1)

        # Evaluate fitness function
        # This is the slowest part, which I'm speeding up by using
        # numexpr. It provides a ~40% gain in execution speed.

        fit_vec = numexpr_evaluate('''N_k * log(N_k/ T_k) ''',
                                   optimization='aggressive')

        p = priors[R]

        A_R = fit_vec - p

        A_R[1:] += best[:R]

        i_max = argmax(A_R)

        last[R] = i_max
        best[R] = A_R[i_max]

        # if(myLikelihood):
        #  logger.debug("Maximum old: %i, Maximum new: %i" %(i_max,i_max_new))
        #  logger.debug("Best old: %s, Best new: %s" %(best[R],best_new[R]))

    pass

    numexpr.set_vml_accuracy_mode(oldaccuracy)

    # if(myLikelihood):
    #   from operator import itemgetter
    #   index, element = max(enumerate(TSs), key=itemgetter(1))
    #   t1,t2 = times[index]
    #   print("Maximum TS is %s in time interval %s-%s" %(element,t1,t2))
    #
    #   best = best_new
    #   last = last_new

    # map(oneLoop,range(N))

    logger.debug("Done\n")

    # Now find blocks
    change_points = np.zeros(N, dtype=int)
    i_cp = N
    ind = N
    while True:
        i_cp -= 1
        change_points[i_cp] = ind

        if ind == 0:
            break

        ind = last[ind - 1]

    change_points = change_points[i_cp:]

    edg = edges[change_points]

    # Transform the found edges back into the original time system
    if (bkgIntegralDistr is not None):
        finalEdges = map(lambda x: lookupTable[x], edg)
    else:
        finalEdges = edg
    pass

    return np.asarray(finalEdges)


# To be run with a profiler
if __name__ == "__main__":

    tt = np.random.uniform(0, 1000, int(sys.argv[1]))
    tt.sort()

    with open("sim.txt", "w+") as f:
        for t in tt:
            f.write("%s\n" % (t))

    res = bayesian_blocks(tt, 0, 1000, 1e-3, None)
    print res
