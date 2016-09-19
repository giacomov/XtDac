#!/usr/bin/env python

# Author: Giacomo Vianello (giacomov@stanford.edu)

# Fix the Good Time Intervals for the XMM-Newton PN instrument

import argparse
import os

try:
  import astropy.io.fits as pyfits
except:
  #If this fail there is no way out
  import pyfits
pass

import numpy
import logging
import warnings
import re
from xtwp4.BayesianBlocks import bayesian_blocks

import multiprocessing

def validFITSfile(arg):
    if not os.path.isfile(arg):
        log.error("The file %s does not exist!" % arg)
    
    #Try to open it to see if it is a FITS file
    try:
      pyfits.open(arg)
    except:
      log.error("The file %s is not a valid FITS file" % arg)
    return arg

def worker(thisBunchEvts):
    boundaries              = []
        
    thisIntervals           = bayesian_blocks(thisBunchEvts,1e-6)
    
    for t1,t2 in zip(thisIntervals[:-1],thisIntervals[1:]):
      
      idx                   = (thisBunchEvts > t1) & (thisBunchEvts < t2)
      
      if(numpy.sum(idx)==2):
        gap_start           = thisBunchEvts[idx][0]+1e-3
        gap_stop            = thisBunchEvts[idx][1]-1e-3
        log.warning("Found a gap ( %s - %s )" %(gap_start,gap_stop))
        boundaries.extend([gap_start,gap_stop])
      pass
    
    pass
    return boundaries


parser                        = argparse.ArgumentParser(description="Fix the GTI accounting for gaps due to counting mode (if any)")

parser.add_argument("-e","--eventfile",help="Event file",
                    type=validFITSfile,required=True)
parser.add_argument("-o","--outfile",help="Name for the output GTI file",required=True)
parser.add_argument("-q","--quadrant",help="Quandrant to use (1,2,3 or 4)",
                                    required=True,type=int,choices=[1,2,3,4])
parser.add_argument("-b","--binsize",help="Binsize for the initial light curve (default: 15)",
                                    required=False,default=15.0,type=float)
parser.add_argument("-m","--mode",help="Mode: 'bblock' for Bayesian block (suggested), or 'fixed' (faster but less accurate)",
                                  required=False,default='bblock',type=str,choices=['bblock','fixed'])
parser.add_argument("-v","--verbosity",required=False,default='info',
                    choices=['info','debug'])

#Main code
if __name__=="__main__":
  
  #Parse arguments
  args                          = parser.parse_args()
  
  #Set up the logger
  levels                        = {'info': logging.INFO, 'debug': logging.DEBUG}
  logging.basicConfig(level=levels[args.verbosity])
  log                           = logging.getLogger("xtfixgti")
  
  ccdmins                       = [1,4,7,10]
  ccdmaxs                       = [3,6,9,12]
  
  ccdmin                        = ccdmins[args.quadrant-1]
  ccdmax                        = ccdmaxs[args.quadrant-1]

  with pyfits.open(args.eventfile,memmap=True) as evt:
      ccdnr                     = evt['EVENTS'].data.field("CCDNR")
      idx                       = ( ccdnr >= ccdmin) & (ccdnr <= ccdmax)
      times                     = evt['EVENTS'].data[idx].field("TIME")
      times.sort()
            
      ontime                    = evt['EVENTS'].header.get("ONTIME")
      tstart_obs                = evt['EVENTS'].header.get("TSTART")
      tstop_obs                 = evt['EVENTS'].header.get("TSTOP")
  pass
  
  log.info("Found %s events in file %s in CCD %i - %i (quadrant %s)" %(times.shape[0],args.eventfile,ccdmin,ccdmax,args.quadrant))
  
  if(args.mode=='bblock'):
    #Read GTIs
    
    #GTIs for different CCDs in this quadrant
    #can be slightly different, read them all and then take as GTI
    #the time interval between the minimum tstart and the maximum
    #tstop
    
    stdgtis                     = []
    for i in range(ccdmin,ccdmax+1):
      stdgtis.append(pyfits.getdata(args.eventfile,"STDGTI%02i" %(i)))
    pass
    
    gtis                        = []
    gtiduration                 = 0
    for i in range(stdgtis[0].field("START").shape[0]):
      tstarts                   = map(lambda x:x.field("START")[i],stdgtis)
      minStart                  = min(tstarts)
      
      tstops                    = map(lambda x:x.field("STOP")[i],stdgtis)
      maxStop                   = max(tstops)
      
      gtis.append((minStart,maxStop))
      gtiduration              += maxStop-minStart
    pass
    
    #Now divide data in bunches and run the Bayesian Block
    #algorithm on each bunch   
    #Note that the "bunches" list will contain references
    #to the data, not copies, so the memory usage will
    #not grow
    
    bunches                     = []
    
    for i,(tstart,tstop) in enumerate(gtis):
      
      #Select events between tstart and tstop
      idx                       = (times >= tstart) & (times <= tstop)
      thisEvts                  = times[idx]
            
      #Divide in bunches of bunchSize events to speed up computation  
      bunchSize                 = 15000 #events
      nBunches                  = int(numpy.ceil(thisEvts.shape[0]/float(bunchSize)))      
      
      newBunches                = []
      for bunch in range(nBunches):
        
        newBunches.append(thisEvts[bunch*bunchSize:min((bunch+1)*bunchSize,thisEvts.shape[0]-1)])
      
      pass
      
      if(len(newBunches)>1):
        #Add also bunches shifted by half period
        for bunch in range(nBunches):
          
          newBunch              = thisEvts[bunch*bunchSize+bunchSize/2:
                                     min((bunch+1)*bunchSize+bunchSize/2,
                                     thisEvts.shape[0]-1)]
          if(len(newBunch)>0):
            newBunches.append(newBunch)
        
        pass
      pass
        
      pass
      
      bunches.extend(newBunches)
      
    pass

    #chunkSize                 = 10000 #Number of events for each worker
    #nDivisions                = int(numpy.ceil(times.shape[0]/chunkSize))
    
    pool                      = multiprocessing.Pool(8)
    results                   = pool.map(func=worker, iterable=bunches)
    pool.close()
    
    #Now process the results and reconstruct the list of
    #Bad Time Intervals (BTIs)
        
    BTIs                      = []
    
    for r in results:
      if(len(r)==0):
        #No bad time interval for this chunk
        continue
      else:
        #Bad time intervals for this chunk
        
        #Finally append the gaps for this chunk to the list of bad time intervals
        #boundaries
        BTIs.extend(r)
      pass
    pass
    
    BTIs                      = numpy.asarray(BTIs)
    BTIs.sort()
    
    #Remove duplicates (possibly produced by the shifted search above)
    BTIs                      = numpy.unique(BTIs)
    
    #Now build the new GTIs, which are the old GTIs minus
    #the time of the gaps found in the previous passages
    
    tstarts                   = []
    tstops                    = []
    gaps                      = 0
    lostTime                  = 0
    
    for (t1,t2) in gtis:
      
      tstarts.append(t1)
      
      #Find the gaps within this GTI
      idx                     = (BTIs > t1) & (BTIs < t2)
      if(numpy.sum(idx)==0):
        continue
      else:
        bti_tstarts           = BTIs[idx][:-1]
        bti_tstops            = BTIs[idx][1:]
        
        assert(len(bti_tstarts)==len(bti_tstops))
        
        tstarts.extend(bti_tstops)
        tstops.extend(bti_tstarts)
        
        #Increment the counter of the number of gaps,
        #compute the duration of the gaps for this chunk and
        #increment the total lost time counter accordingly
        gaps                 += len(bti_tstops)
        durations             = map(lambda x:x[1]-x[0],zip(bti_tstarts,bti_tstops))
        lostTime             += numpy.sum(durations)
      pass
      
      tstops.append(t2)
    pass
            
    log.info("Number of GTIs already present: %i" %(len(gtis)))
    log.info("Total time in GTIs before analysis: %.2f" %(gtiduration))
    log.info("Found %i gaps not accounted for in GTIs, for a total of %.2f seconds of dead time" %(gaps,lostTime))
    log.info("New total time in GTIs: %.2f" %(gtiduration-lostTime))
        
  elif(args.mode=='fixed'):
  
    #Make the initial light curve with fixed bin
    binsize                     = args.binsize
    Nbins                       = int(numpy.ceil((tstop_obs-tstart_obs)/binsize))
    res                         = numpy.histogram(times,Nbins)
    
    centers                     = (res[1][1:]+res[1][:-1])/2.0
    rates                       = res[0]/binsize
    
    idx                         = (rates==0)
    centersOfEmptyBins          = centers[idx]
    
    #For each bin with zero counts optimize the boundaries
    newBTIs = []
    
    for c in centersOfEmptyBins:
        
        lowerBound              = c - (binsize/2.0)
        
        #Find the closest event before lowerBound (look at numpy.searchsorted manual
        #to undertstand this)
        evid                    = numpy.searchsorted(times,lowerBound)-1
        
        evtBefore               = times[evid]
        log.debug("Found event at time %s before boundary %s" %(evtBefore,lowerBound))
        
        upperBound              = c + binsize/2.0
        
        #Find the closest event after upperBound
        evid                    = numpy.searchsorted(times,upperBound)
        
        #If the observation ends with a Bad Time Interval, there are no
        #events after the upperBound, thus evid will be larger than
        #the last element in the array
        if(evid==times.shape[0]):
          #Use the stop time of the observation as end of this BTI
          evtAfter              = tstop_obs
        else:
          evtAfter              = times[evid]
        pass
        
        log.debug("Found event at time %s after boundary %s" %(evtAfter,upperBound))
        
        #Append new Bad Time Interval
        newBTIs.append((evtBefore+1e-2,evtAfter-1e-2))
        log.debug("Added BTI %s - %s" %(newBTIs[-1][0],newBTIs[-1][1]))
    pass
    
    #Now collapse the BTIs
    newBTIs                     = set(newBTIs)
    
    log.info("Found %s Bad Time Intervals" %(len(newBTIs)))
  
    
    #Now generate the GTIs
    #The first GTI starts at TSTART (beginning of observation),
    #but we will add it later
    tstarts                     = [tstart_obs]
    #All the other start time of GTIs are the stop time of the
    #Bad Time Intervals
    tstarts.extend(map(lambda x:x[1],newBTIs))
    #Convert to numpy (needed for pyfits)
    tstarts                     = numpy.asarray(tstarts)
    
    #All the stop times of the GTIs are the start time of the
    #Bad Time Intervals
    tstops                      = map(lambda x:x[0],newBTIs)
    
    #Last tstop is the time of the last event
    tstops.append(times.max()+1e-2)
    
    #Convert to numpy (needed for pyfits)
    tstops                      = numpy.asarray(tstops)
    
    idx                         = numpy.argsort(tstarts)
    tstarts                     = tstarts[idx]
    tstops                      = tstops[idx]
  pass
    
  #generate the FITS file with the GTIs
  c1                          = pyfits.Column(name='START', format='D', unit='S', array=tstarts)
  c2                          = pyfits.Column(name='STOP', format='D', unit='S', array=tstops)
  coldefs                     = pyfits.ColDefs([c1, c2])
  tbhdu                       = pyfits.BinTableHDU.from_columns(coldefs)
  
  #Prepare the header to conform to XMM standards
  header                      = tbhdu.header
  header.set("EXTNAME","STDGTI")
  header.set("HDUCLASS","OGIP")
  header.set("HDUCLAS1","GTI")
  header.set("HDUCLAS2","STANDARD")
  header.set("ONTIME",ontime)
  header.set("TSTART",tstart_obs)
  header.set("TSTOP",tstop_obs)
  header.set("TIMEUNIT",'s')
  header.set("TIMESYS","TT")
  header.set("MJDREF",50814.0)
  header.set("TIMEREF","LOCAL")
  header.set("TASSIGN","SATELLITE")
  header.set("TIMEZERO",0)
  header.set("CLOCKAPP",True)
  header.set("MINGTISZ",0)
  header.set("SUMSHORT",0)
  
  #Finally write the file
  try:
    tbhdu.writeto(args.outfile,clobber=True)
  except:
    log.fatal("Could not write outfile. Disk full? Problems with write permission in %s ?" %(os.getcwd()))
    raise
  else:
    log.info("%i GTIs written in %s" %(len(tstarts),args.outfile))
pass
