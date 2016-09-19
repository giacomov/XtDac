# Author: Giacomo Vianello (giacomov@stanford.edu)

from scipy import interpolate
from sklearn.neighbors import KernelDensity
from astroML.density_estimation import KNeighborsDensity
import numpy
import matplotlib.pyplot as plt
import sys

class DensityEstimation(object):
    def __init__(self):
        raise NotImplemented("This class must be subclassed.")
        
    def plot(self,binsize):
    
        fig                        = plt.figure()
        sub                        = fig.add_subplot(111)
        
        #Data histogram
        nbins                      = int((self.dt)/binsize)
        n,bins,patches             = sub.hist(self.evts, nbins,alpha=0.5,normed=False)
        
        #Retrieve the effective binsize
        binsize                    = bins[1]-bins[0]
        
        #Model plot
        x_grid                     = numpy.arange(self.t1,self.t2,10.0)
        sub.plot(x_grid, self.model(x_grid)*self.evts.shape[0]*binsize, '--',color='blue', alpha=1, lw=1,zorder=10)
        sub.set_xlim([self.t1,self.t2])
        
        return fig,sub
    
    def getIntegralDistribution(self,t0=0,renorm=1):
        
        #Define the x grid
        gridsize                   = min(5,(self.dt)/2)
        x_grid                     = numpy.arange(self.t1-t0,self.t2-t0,gridsize)
        
        #Add a point before t1 with 1e-5 counts, so that extrapolating
        #backwards will return 1e-5, instead of negative or zero counts
        x_grid                     = numpy.insert(x_grid,0,self.t1-t0-1e-3)
        
        #Compute the integral distribution for each point in the x grid
        y                          = [0]
        for x in x_grid[1:]:
            y.append(self.model.integral(self.t1,x+t0)*self.evts.shape[0]*renorm)
        
        #Build a linear spline interpolating the integral distribution,
        #which is very fast to evaluate on any given x
        intDistr                   = interpolate.InterpolatedUnivariateSpline(x_grid,y,k=1)
        
        return intDistr


class KDE(DensityEstimation):
    def __init__(self,evts,t1,t2,bw=50):    
        
        self.evts                  = evts
        
        #Decide the bandwidth
        self.t1                    = t1
        self.t2                    = t2
        self.dt                    = self.t2-self.t1
        bw                         = min(bw,int(self.dt))
        
        #Perform the Kernel Density Estimation
        self.bandwidth             = bw
        kde                        = KernelDensity(bandwidth=bw, kernel='gaussian',rtol=0.001)
        #sys.stderr.write("Fitting...")
        kde.fit(evts[:, numpy.newaxis])
        #sys.stderr.write("done")
        
        #Evaluate the KDE and interpolate it
        #sys.stderr.write("Evaluating...")
        x_grid                     = numpy.arange(self.t1,self.t2,10.0)
        pdf                        = numpy.exp(kde.score_samples(x_grid[:, numpy.newaxis]))
        #sys.stderr.write("done")
        
        #sys.stderr.write("Interpolating...")
        self.model                 = interpolate.InterpolatedUnivariateSpline(x_grid,pdf,k=1)
        #sys.stderr.write("done")
        
pass

class KNeighbors(DensityEstimation):
   def __init__(self,evts,t1,t2,k=10000):
        self.evts                  = evts
        
        #Decide the bandwidth
        self.t1                    = t1
        self.t2                    = t2
        self.dt                    = self.t2-self.t1
        
        #Perform the Kernel Density Estimation
        knd                        = KNeighborsDensity('bayesian', n_neighbors=k)
        sys.stderr.write("Fitting...")
        knd.fit(evts[:, numpy.newaxis])
        sys.stderr.write("done")
        
        #Evaluate the KDE and interpolate it
        sys.stderr.write("Evaluating...")
        x_grid                     = numpy.arange(self.t1,self.t2,10.0)
        pdf                        = knd.eval(x_grid[:, numpy.newaxis])/evts.shape[0]
        sys.stderr.write("done")
        
        sys.stderr.write("Interpolating...")
        self.model                 = interpolate.InterpolatedUnivariateSpline(x_grid,pdf,k=1)
        sys.stderr.write("done")
        
