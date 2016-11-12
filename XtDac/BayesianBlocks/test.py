#!/usr/bin/env python

from BayesianBlocks import *
import numpy as np
import sys

# To be run with a profiler
if __name__ == "__main__":
    
    tt = np.random.uniform(0, 1000, int( sys.argv[1] ))
    tt.sort()

    with open("sim.txt","w+") as f:
        for t in tt:
            f.write("%s\n" %(t))

    res = bayesian_blocks(tt, 0, 1000, 1e-3, None)
    print res
