from astropy.io import fits as pyfits

import pyregion

import warnings

import numpy

class FitsRegionFile(object):
    
    def __init__(self, filename, minimumSize=0):
        
        with pyfits.open(filename) as f:
            
            header = f['EVENTS'].header
            
            self.X = f['EVENTS'].data.field("X")
            self.Y = f['EVENTS'].data.field("Y")
            self.t = f['EVENTS'].data.field("TIME")
            
            self.regions = []
            
            for i,row in enumerate( f['REGION'].data ):
                
                X = row.field("X")[0]
                Y = row.field("Y")[0]
                r1,r2 = row.field("R")[:2]
                rot = row.field("ROTANG")[0]
                
                reg = "physical;ellipse({X},{Y},{r1},{r2},{rot})".format(**locals())
                
                try:
                    
                     r = pyregion.parse(reg).as_imagecoord( header )
                
                except:
                    
                    warnings.warn("Could not parse region %s" %(i+1))
                
                if r1 >= minimumSize and r2 >= minimumSize:
                    
                    self.regions.append( r )
                
                else:
                    
                    warnings.warn("Removing region %s because is too small" %(i+1))
                    
    
    def iteritems(self):
        
        '''Loop over the non-empty regions returning the corresponding events. Use like this:
           
           rr = FitsRegionFile( filename )
           
           for x,y,t,filt,reg in rr.iteritems():
               
               ...
        
        '''
                
        for i in range( len( self.regions ) ):
            
            filt = self.regions[i].get_filter()
            
            res = filt.inside( self.X, self.Y )
            
            if self.X[res].shape[0] == 0:
                
                #Empty. Skip to the next one
                
                continue
            
            yield self.X[res], self.Y[res], self.t[res], filt, self.regions[i]
            

class FitsRegion( FitsRegionFile ):
    
    def __init__(self, X, Y, t, header, regionDefinition):
        
        self.X = X
        self.Y = Y
        self.t = t
        
        self.regions = [ pyregion.parse( regionDefinition ).as_imagecoord( header ) ]
