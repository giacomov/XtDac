#!/usr/bin/env python

#Author: Giacomo Vianello (giacomov@stanford.edu)

##############
# Produce background maps for XMM EPIC via smart smoothing
##############


import argparse
import os, sys
import numpy
import logging
import warnings

import scipy.fftpack

import astropy.io.fits as pyfits
from astropy.convolution import Gaussian2DKernel

def validFITSfile(arg):
    if not os.path.isfile(arg):
        log.error("The file %s does not exist!" % arg)
    
    #Try to open it to see if it is a FITS file
    try:
      pyfits.open(arg)
    except:
      log.error("The file %s is not a valid FITS file" % arg)
    return arg

def smooth(imgtot,sigma=100):
    kernel = Gaussian2DKernel(sigma)
    kernel.normalize()
        
    result = convolve_fft(imgtot, kernel,interpolate_nan=True,
                          ignore_edge_zeros=True, 
                          complex_dtype='complex64')
    return result


parser                        = argparse.ArgumentParser(description="Create a background map for an XMM image")

parser.add_argument("-i","--image",help="Cleaned image",type=validFITSfile,required=True)
parser.add_argument("-e","--expomap",help="Exposure map (do not remove point source regions)",
                    type=validFITSfile,required=True)
parser.add_argument("-m","--mask", help="Cheese mask",type=validFITSfile,required=True)
parser.add_argument("-o","--outfile",help="Name for the output Background map",required=True)
parser.add_argument("-s","--smoothing", 
                    help="Sigma of the Gaussian kernel used for the smoothing (default: 80)",
                    type=int, default=80,required=False)
parser.add_argument("-v","--verbosity",required=False,default='info',
                    choices=['info','debug'])

def go(args):  
  #Set up the logger
  levels                        = {'info': logging.INFO, 'debug': logging.DEBUG}
  logging.basicConfig(level=levels[args.verbosity])
  log                           = logging.getLogger("xtbgmap")

  
  #Read in images
  expot                         = pyfits.getdata(args.expomap,0)
  log.debug("Read an exposure map with shape %s and %s non-null elements" %(str(expot.shape),numpy.sum((expot!=0))))
  
  maskt                         = pyfits.getdata(args.mask,0)
  log.debug("Read a mask with shape %s and %s non-null elements" %(str(maskt.shape),numpy.sum((maskt!=0))))
  
  #Cast the image to float so we can use numpy.nan (and divide by the exposure
  #without worrying about round errors)
  img                           = numpy.array(pyfits.getdata(args.image,0),'f')
  #Get header to keep the WCS of the original image
  head                          = pyfits.getheader(args.image,0)
  log.debug("Read an image with shape %s and %s non-null elements" %(str(img.shape),numpy.sum((img!=0))))
    
  #Makes all the pixels with exposure=0 as NaN in the image
  #so that they will filled with interpolated values before
  #the smoothing
  idx                           = (expot==0)
  img[idx]                      = numpy.nan
  log.debug("Changed %s pixels with zero exposure to NaN in the image" %(numpy.sum(idx)))
  
  #Now mask out the pixels which are empty in the mask
  #(but are not zero in the exposure)
  idx2                          = (maskt==0)
  img[idx2]                     = numpy.nan
  log.debug("Changed %s pixels to NaN in the image to reflect the mask" %(numpy.sum(idx2)))
  
  log.info("Smoothing...")
  with numpy.errstate(all='ignore'):
    newImg                        = smooth(img,args.smoothing)
  log.info("done")
    
  #Correct for the exposure
  log.debug("Dividing the smoothed image by the exposure...")
  idx3                          = (expot > 0)
  newImg[idx3]                 /= expot[idx3]
  log.debug("done")
  
  #Massage the final output
  
  #Restore all values with exposure 0 as zeroes (instead of NaN)
  log.debug("Changing back to zero all pixels with zero exposure")
  newImg[idx]                   = 0
  
  #Write the output file
  log.debug("Writing output file...")
  try:
    pyfits.writeto(args.outfile,newImg,header=head,clobber=True)
  except:
    raise RuntimeError("Could not write output file. Check free space or write permissions for output directory.")
  else:
    log.debug("done")
  pass

#The following is taken from astropy.convolution.convolve_fft,
#but has been optimized for this particular case, removing some unused options
#and optimizing memory usage. The original function was using several Gb of memory
#This implementation uses >4x less memory
def convolve_fft(array, kernel, fill_value=0, crop=True,
                 return_fft=False, fft_pad=True, psf_pad=False,
                 interpolate_nan=False, quiet=False, ignore_edge_zeros=False,
                 min_wt=0.0, normalize_kernel=False, allow_huge=False,
                 fftn=numpy.fft.fftn, ifftn=numpy.fft.ifftn,
                 complex_dtype=numpy.complex):

    kernel = kernel.array

    if numpy.abs(kernel.sum() - 1) < 1e-8:
        kernel_is_normalized = True
    else:
        raise RuntimeError("Kernel must be normalized")
   
    # Check that the number of dimensions is compatible
    if array.ndim != kernel.ndim:
        raise ValueError("Image and kernel must have same number of "
                         "dimensions")

    arrayshape = array.shape
    kernshape = kernel.shape
    
    # NaN and inf catching
    nanmaskarray = numpy.isnan(array) + numpy.isinf(array)
    array[nanmaskarray] = 0
    
    fsize = 2 ** numpy.ceil(numpy.log2(
                   numpy.max(numpy.array(arrayshape) + numpy.array(kernshape))))
    newshape = numpy.array([fsize for ii in range(array.ndim)], dtype=int)
    
    # separate each dimension by the padding size...  this is to determine the
    # appropriate slice size to get back to the input dimensions
    arrayslices = []
    kernslices = []
    for ii, (newdimsize, arraydimsize, kerndimsize) in enumerate(
                                         zip(newshape, arrayshape, kernshape)):
        center = newdimsize - (newdimsize + 1) // 2
        arrayslices += [slice(center - arraydimsize // 2,
                              center + (arraydimsize + 1) // 2)]
        kernslices += [slice(center - kerndimsize // 2,
                             center + (kerndimsize + 1) // 2)]

    if not numpy.all(newshape == arrayshape):
        bigarray = numpy.ones(newshape, dtype=complex_dtype) * fill_value
        bigarray[arrayslices] = array
    else:
        bigarray = array

    if not numpy.all(newshape == kernshape):
        bigkernel = numpy.zeros(newshape, dtype=complex_dtype)
        bigkernel[kernslices] = kernel
    else:
        bigkernel = kernel
    
    arrayfft = fftn(bigarray)
    # need to shift the kernel so that, e.g., [0,0,1,0] -> [1,0,0,0] = unity
    kernfft = fftn(numpy.fft.ifftshift(bigkernel))
    fftmult = arrayfft * kernfft
    
    #Explicitly clean memory to avoid high peak-memory usage
    del arrayfft
    del bigkernel
    del bigarray
        
    if (interpolate_nan or ignore_edge_zeros) and kernel_is_normalized:
        if ignore_edge_zeros:
            bigimwt = numpy.zeros(newshape, dtype=complex_dtype)
        else:
            bigimwt = numpy.ones(newshape, dtype=complex_dtype)
        bigimwt[arrayslices] = 1.0 - nanmaskarray * interpolate_nan
        wtfft = fftn(bigimwt)
        # I think this one HAS to be normalized (i.e., the weights can't be
        # computed with a non-normalized kernel)
        wtfftmult = wtfft * kernfft / kernel.sum()
        
        del wtfft
        del kernfft
        
        wtsm = ifftn(wtfftmult)
        
        del wtfftmult
        
        # need to re-zero weights outside of the image (if it is padded, we
        # still don't weight those regions)
        bigimwt[arrayslices] = wtsm.real[arrayslices]
        
        del wtsm
        
        # curiously, at the floating-point limit, can get slightly negative numbers
        # they break the min_wt=0 "flag" and must therefore be removed
        bigimwt[bigimwt < 0] = 0
    else:
        bigimwt = 1

    if numpy.isnan(fftmult).any():
        # this check should be unnecessary; call it an insanity check
        raise ValueError("Encountered NaNs in convolve.  This is disallowed.")

    # restore NaNs in original image (they were modified inplace earlier)
    # We don't have to worry about masked arrays - if input was masked, it was
    # copied
    array[nanmaskarray] = numpy.nan
    #kernel[nanmaskkernel] = numpy.nan

    if return_fft:
        return fftmult

    if interpolate_nan or ignore_edge_zeros:
        rifft = (ifftn(fftmult)) / bigimwt
        if not numpy.isscalar(bigimwt):
            rifft[bigimwt < min_wt] = numpy.nan
            if min_wt == 0.0:
                rifft[bigimwt == 0.0] = 0.0
    else:
        rifft = (ifftn(fftmult))

    if crop:
        result = rifft[arrayslices].real
        return result
    else:
        return rifft.real


#Main code
if __name__=="__main__":
  #Parse arguments
  args                          = parser.parse_args()
  go(args)
