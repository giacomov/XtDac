import astropy.io.fits as pyfits
from astropy.wcs import WCS
import numpy
import sys
import scipy.interpolate

class Exposure(object):
  
  def __init__(self,exposureMap):
    
    with pyfits.open(exposureMap) as f:
      self.exposureMap_ = f[0].data
      
      self.header = f[0].header
        
    self.exposureWCS = WCS(self.header)
    
    #Interpolate exposure map
    self.exposureMap = scipy.interpolate.RectBivariateSpline(
                           range(self.exposureMap_.shape[0]),
                           range(self.exposureMap_.shape[1]),
                           self.exposureMap_,kx=1,ky=1,s=0)
    
    self.onTime = self.header['EXPOSURE']
  
  def getMaximumExposure(self):
    
    return self.exposureMap_.max()
  
  def getAverageExposure(self):
    return self.averageExposure
  
  def getExposureAtSky(self,ras,decs):
    i,j = self.exposureWCS.wcs_world2pix(ras,decs,0)
        
    return self.exposureMap.ev(j,i)
  
  def cacheShape(self,xmin,xmax,ymin,ymax,header):
    #This must be called before attempting to use the
    #getExposureAtXY method
    
    #Build WCS from input header
    guestWCS = WCS(naxis=2)
    guestWCS.wcs.crpix =[header['REFXCRPX'],header['REFYCRPX']]    
    guestWCS.wcs.cdelt = [header['REFXCDLT'],header['REFYCDLT']]
    guestWCS.wcs.crval = [header['REFXCDLT'],header['REFYCDLT']]   
    guestWCS.wcs.crval = [header['REFXCRVL'],header['REFYCRVL']]
    guestWCS.wcs.ctype = [header['REFXCTYP'],header['REFYCTYP']]
    
    self.guestWCS = guestWCS
    
    #Get the minimum and maximum X and Y contained
    #in the image/event file
    xxmin = header['TDMIN6']
    xxmax = header['TDMAX6']
    yymin = header['TDMIN7']
    yymax = header['TDMAX7']
    
    if(xmin < xxmin or xmax > xxmax or
       ymin < yymin or ymax > yymax):
       raise RuntimeError("Requested min/max value outside of image")
        
    #cache the exposure for each coordinate in x and y
    x = numpy.arange(xmin,xmax,1)
    y = numpy.arange(ymin,ymax,1)
    
    xx,yy = numpy.meshgrid(x,y,indexing='ij')
    
    ras,decs = guestWCS.wcs_pix2world(xx,yy,1)
        
    self.exposures_ = self.getExposureAtSky(ras,
                                           decs)
    
    idx = self.exposures_ <= 0
    
    self.averageExposure = numpy.median(self.exposures_[~idx].flatten())
    
    #self.exposures_[idx] = self.averageExposure / 1e3
    
    self.exposures = scipy.interpolate.RectBivariateSpline(
                                       x,y,self.exposures_,
                                       kx=1,ky=1,s=0)
        
    

    
  def getOnTime(self):
    return self.onTime
  
  def getExposureAtXY(self,x,y):
    return self.exposures.ev(x,y)
