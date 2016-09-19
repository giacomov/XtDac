#!/usr/bin/env python

# Author: Giacomo Vianello (giacomov@stanford.edu)

import argparse
from xtwp4.xmmdac import TransientDiscovery
from xtwp4.xmmdac import AntiCoincidence
import os

#Parse input arguments

parser                        = argparse.ArgumentParser()

parser.add_argument("eventFile",help="Event file to use (already cleaned and selected in energy)")

parser.add_argument("outfile",help="Outfile which will contain the time intervals")

parser.add_argument("-m","--maxDuration",help="Maximum duration of intervals to keep (default=1000 s)", 
                    type=float,default=1000)

parser.add_argument("-c","--ncpus", help="Number of CPUs to use (default=1)", type=int, default=1)

parser.add_argument("-t","--acdTolerance",help="Tolerance in seconds for the anti-coincidence operation."+ 
                                               " If two intervals starts within this tolerance, and they are"+
                                               " not from two close regions, they are considered background"+
                                               " variations and are excluded from the results", type=float,
                                               default=10) 

parser.add_argument("-f","--falsePositive",help="False positive rate for the Bayesian Blocks algorithm."+ 
                                               " Note that this is pre-trial, and for a single cell.", type=float,
                                               default=0.01) 


parser.add_argument("-v","--verbose",help="Print the time intervals found by the analysis",
                    action="store_true")

args                          = parser.parse_args()

#Main code
if __name__=="__main__":
  print("Using %s CPUs in processing event file %s" %(args.ncpus,args.eventFile))
  t                           = TransientDiscovery.TransientDiscovery(args.eventFile)
  
  print("Dividing events in regions:")
  t.preSelect()
  print("done")
  
  print("Looking for excesses...")
  intervalsBunches            = t.goMulti(args.falsePositive,args.ncpus)
  print("done")
  
  with open("__intervals.txt","w+") as f:
    f.write("#Tstart Tstop Duration Box_corner_X Box_corner_Y Box_width Box_height Rate Excess\n")
    for bunch in intervalsBunches:
      if(bunch==None):
        continue
      for interv in bunch:
        f.write("%s %s %s %s %s %s %s %s %s\n" %(interv.tstart,interv.tstop,interv.tstop-interv.tstart,
                                        interv.region.c1,interv.region.c2,
                                        interv.region.width,interv.region.height,interv.rate,interv.excess))
      pass
    pass
  pass
  
  acd                         = AntiCoincidence.AntiCoincidence("__intervals.txt")
  intervals                   = acd.go(args.acdTolerance)
  #os.remove("__intervals.txt")
  
  idx                         = (intervals[2] < args.maxDuration)
  print("Found %s intervals smaller than %s seconds\n" %(intervals[0][idx].shape[0],args.maxDuration))
  with open(args.outfile,'w+') as f:
    for t1,t2 in zip(intervals[0][idx],intervals[1][idx]):
      f.write("%s %s\n" %(t1,t2))
      if(args.verbose):
        print("%s %s" %(t1,t2))
    pass
  pass
  print("written in %s" %(args.outfile))
