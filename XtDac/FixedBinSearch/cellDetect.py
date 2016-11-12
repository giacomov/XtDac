import pyfits

import scipy.stats

import numpy

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import math

import numexpr

import cartesian

import logging

log = logging.getLogger("main")

class Cell(object):
    
    def __init__(self, x,y, w, h):
        
        #x,y are the coordinates of the upper left corner
        self.x = x
        self.y = y
        self.w = w
        self.h = h
        
        self.cx = self.x + self.w/2.0
        self.cy = self.y + self.h/2.0
    
    def __repr__(self):
        
        rep = ( "(x,y) = (%s, %s), w = %s, h = %s, (xc,yc) = (%s, %s)" 
                % (self.x, self.y, self.w, self.h, self.cx, self.cy) )
        
        return rep
    
    def count(self, X,Y):
        idx = ( (X >= self.x) & (X <= self.x + self.w) 
              & (Y >= self.y) & (Y <= self.y + self.h) )
        self.counts = numpy.sum( idx   )
        
        return idx
    
    def getCenter(self):
        
        return [ self.cx, self.cy ]
    
    def setSignificance(self, B, s):
        
        #Significance in presence of background uncertainty
        
        o = self.counts
        
        if o < B:
           
            s = -1
        
        else:
            
            b0 = 0.5 * ( math.sqrt( pow(B,2) - 2 * s * (B - 2 * o) + pow(s,4) ) + B - pow(s,2) )
            
            s = math.sqrt(2) * math.sqrt( o * math.log( o / b0 ) + ( pow(b0 - B, 2) / (2 * s) ) + b0 - o )
        
        self.significance = s
    
    def setProbability(self,p):
        
        self.prob = p
    
    def distance(self, x,y):
        
        return math.sqrt( pow(x-self.cx, 2) + pow(y-self.cy, 2) )
    
class CellDetect(object):
    
    def __init__(self, X, xmin, xmax, Y, ymin, ymax, cellsize, multiplicity=2):
        
        self.X = numpy.asarray( X )
        self.xmin = xmin
        self.xmax = xmax
        
        self.Y = numpy.asarray( Y )
        self.ymin = ymin
        self.ymax = ymax
        
        assert self.X.shape[0] == self.Y.shape[0]
        
        self.cellsize = cellsize
        
        self.gridGen( multiplicity )
        
        log.debug("CellDetect initiated.")
        log.debug("Number of events: %s" %(self.X.shape[0]))
        
    def gridGen(self, multiplicity ):
        
        self.multiplicity = int( multiplicity )
        
        maxx = self.xmax
        minx = self.xmin
        
        maxy = self.ymax
        miny = self.ymin
        
        approxSize = self.cellsize
        
        #Adjust the horizonthal size so that we have an integer number
        #of boxes
        width                     = (maxx-minx)
        #The division and multiplication by 2 is a trick to make sure that
        #the width divided by horizSize/2 is an integer 
        sw                        = numpy.rint(width/approxSize/2.0)*2.0
        horizSize                 = width/sw
    
        #Adjust the vertical size so that we have an integer number
        #of boxes
        height                    = (maxy-miny)
        sw                        = numpy.rint(height/approxSize/2.0)*2.0
        vertSize                  = height/sw
        #Generate boxes
        xcoords                   = numpy.arange(minx,maxx,horizSize/multiplicity)
        ycoords                   = numpy.arange(miny,maxy,vertSize/multiplicity)
        
        cells = []
        
        #Copy the arrays to speed up the filling
        X_ = numpy.array(self.X, copy=True)
        Y_ = numpy.array(self.Y, copy=True)
        
        for x in xcoords:
            
            thisLine = []
            
            for y in ycoords:
                
                thisCell = Cell(x,y,horizSize, vertSize)
                idx = thisCell.count(X_,Y_)
                #X_ = X_[~idx]
                #Y_ = Y_[~idx]
                thisLine.append( thisCell  )
            
            cells.append(thisLine)
        
        self.cells = numpy.array(cells)
    
    def _left(self,j):
        
        return max(0,j-1)
    
    def _right(self,j):
        
        return min(j+1, self.cells.shape[1] - 1 )
    
    def _up(self,i):
        
        return max(0, i-1)
    
    def _down(self,i):
        
        return min(i+1, self.cells.shape[0] - 1 )        
    
    def findExcesses(self, bkgLevel = None, bkgError = None, verbose=False):
        
        cts = numpy.array( map(lambda cell:cell.counts, self.cells.flatten()) )
        
        if bkgLevel is None:
            
            #Auto-determine bkg level
            
            idx = (cts == 0)
            
            #There is a 5% probability of obtaining 0 counts from an expected
            #value of 3. Hence 3 is the 95% upper limit for the count flux
            #if 0 is observed
            cts[idx] = 3.0 / numpy.sum(idx)
            
            if verbose:
              print("Average value: %s" % numpy.average(cts))
              print("Median value: %s" % numpy.median(cts))
              print("Peak value: %s" % cts.max() )
              print("Std. dev: %s" % numpy.std(cts))
            
            bkgLevel = numpy.median(cts)
            bkgError = numpy.std(cts)
        
        else:
            idx = (cts != 0)
            log.debug("Median counts: %s" %(numpy.median(cts[idx])))
            log.debug("Background level: %s" %(bkgLevel))
        
        #pois = scipy.stats.poisson(bkgLevel)
        
        #map(lambda cell:cell.setProbability( pois.sf(cell.counts)+pois.pmf(cell.counts)), self.cells.flatten() )
        
        map(lambda cell:cell.setSignificance( bkgLevel, bkgError ), self.cells.flatten() )
        
        #prob = numpy.zeros_like(self.cells)
        s = numpy.zeros_like( self.cells )
        
        for i,row in enumerate( self.cells ):
            for j,cell in enumerate( row ):
                #prob[i,j] = cell.prob
                s[i,j] = cell.significance
        
        #idx = (prob <= 1e-6)
        
        idx = (s >= 5) #sigma
        
        nexc = numpy.sum(idx)
        
        log.debug("Excesses: %s" % nexc )
        
        if nexc == 0:
           
           self.centroids = numpy.array([],dtype=float)
           
           return self.centroids
        
        ex,ey = idx.nonzero()
        
        sources = []
                
        for i,j in zip(ex,ey):
                #Check if the surrounding cells have a smaller
                #significance
                
                thisP = self.cells[i][j].counts
                
                one = self.cells[i][self._left(j)].counts
                
                two = self.cells[i][self._right(j)].counts
                
                three = self.cells[self._up(i)][j].counts
                
                four = self.cells[self._down(i)][j].counts
                
                five = self.cells[self._up(i)][self._left(j)].counts
                
                six = self.cells[self._up(i)][self._right(j)].counts
                
                seven = self.cells[self._down(i)][self._left(j)].counts
                
                eight = self.cells[self._down(i)][self._right(j)].counts
                
                if ( (one <= thisP) & (two <= thisP) & 
                     (three <= thisP) & (four <= thisP) &
                     (five <= thisP) & (six <= thisP) &
                     (seven <= thisP) & (eight <= thisP)
                   ):
                    #Source!
                    sources.append( self.cells[i][j] )
                    log.debug("Found source with significance %s" %( self.cells[i][j].significance ))
                    
                    if verbose:
                        print([thisP, one,two,three,four,five,six,seven,eight])
        
        if verbose:
            print("")
                
        #Check for duplicated sources (sources in adiacent boxes)
        
        cleanCells = numpy.array( sources, copy=True )
        
        for src in sources:
            
            distances = numpy.zeros( cleanCells.shape[0] )
            
            for i,src2 in enumerate( cleanCells ):
                
                distances[i] = src.distance( *src2.getCenter() )
            
            #Note: the following will select at least 1
            #cell, which is the current one, with distance = 0
            
            idx = ( distances < self.cellsize * math.sqrt(2) )
            
            if distances[idx].shape[0] > 1:
                
                #Duplicate
                log.debug("Warning: confused source. Using barycenter")
                
                #Getting all the centers of the confused cells
                
                xs = []
                ys = []
                
                for src3 in cleanCells[idx]:
                                      
                   x,y = ( src3.x, src3.y )
                   xs.append( x )
                   ys.append( y )
                
                #Computing the average centroid
                
                baryX = numpy.average( numpy.unique( xs ) )
                baryY = numpy.average( numpy.unique( ys ) )
                
                #Removing all the involved cells and add a new one
                #centered in the average centroid
                
                cleanCells = cleanCells[~idx]
                
                newCell = Cell( baryX, baryY, src3.w, src3.h )
                
                cleanCells = numpy.append( cleanCells, newCell )
                
        
        #_ = plt.hist2d(self.X, self.Y, [100,100])
        
        centroids = map( lambda cell: cell.getCenter(), cleanCells )
        
        self.centroids = numpy.array(centroids)
        
        #plt.plot(centroids[:,0],centroids[:,1],'o',markersize=5,color='red')
        
        return self.centroids
    
    def plot(self, xbins, ybins, **kwargs):
        
        fig,sub = plt.subplots(1,1)
        
        _ = sub.hist2d(self.X, self.Y, [xbins, ybins])
        
        if len(self.centroids) > 0:
            
            cc = numpy.array(self.centroids)
            
            sub.scatter(cc[:,0], cc[:,1], s=80,alpha=0.8,
                    facecolors='none', edgecolors='r', lw=2)
        
        if 'filter' in kwargs.keys():
            
            xcs = (xbins[:-1] + xbins[1:] ) / 2.0
            ycs = (ybins[:-1] + ybins[1:] ) / 2.0
            
            points = cartesian.cartesian([xcs, ycs])
            
            idx = kwargs['filter'].inside( points[:,0], points[:,1] )
            
            sub.plot(points[~idx,0], points[~idx,1], 'x', markersize=8, color='white') 
        
        #for cell in self.cells.flatten():
            
        #    sub.add_patch( Rectangle((cell.x, cell.y), cell.w, cell.h, color='red',alpha=0.5, fill=False) )
        
        return fig
