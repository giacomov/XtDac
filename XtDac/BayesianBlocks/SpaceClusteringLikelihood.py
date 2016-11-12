import numpy as np
from iminuit import Minuit
import math
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal
from scipy.special import erf
from scipy.integrate import dblquad
import random

erf = math.erf

def getPoissonVariate(N,tmin,tmax):
    rate = N / (float(tmax)-float(tmin))
    tt = [tmin + random.expovariate(rate)]
    while(True):
        if(tt[-1] > tmax):
            break
        tt.append(tt[-1] + random.expovariate(rate))
        
    tt = np.asarray(tt)
    tt = tt[tt <= tmax]
    
    return tt

def getVoronioCells(t,tmin,tmax):
    edges                          = np.concatenate([[tmin],
                                                    0.5 * (t[1:] + t[:-1]),
                                                    [tmax]])
    return edges

def getIsoLikelihoodValue(N,tmin,tmax):
    return (N * (math.log(N)-math.log(tmax-tmin)))

class UnbinnedModel(object):
    
    def __init__(self, t,tmin,tmax, x, xmin, xmax, y,ymin,ymax):        
        self.t = t       
        self.x = x
        self.y = y
        
        self.N = t.shape[0]
        
        self.setRange(tmin,tmax,xmin,xmax,ymin,ymax,0,self.t.shape[0])
        
        
    def setRange(self,tmin,tmax,xmin,xmax,ymin,ymax,lo_idx,hi_idx):
        self.tmin = tmin
        self.tmax = tmax
        self.dt = tmax-tmin
        
        self.N = hi_idx - lo_idx
        
        self.xmin = xmin
        self.xmax = xmax
        self.dx = self.xmax - self.xmin
        self.NlogDx = self.N * math.log(self.dx)
        
        self.ymin = ymin
        self.ymax = ymax
        self.dy = self.ymax - self.ymin
        self.NlogDy = self.N * math.log(self.dy)
        
        self.dOmega = self.dx * self.dy * self.dt
        
        self.curT = self.t[lo_idx:hi_idx]
        self.curX = self.x[lo_idx:hi_idx]
        self.curY = self.y[lo_idx:hi_idx]
        self.curN = self.curT.shape[0]
                
    def mloglike(self, *parameters):
        
        self.setParameters(*parameters)
        
        firstTerm = np.sum(np.log(self.evaluate(self.curT,self.curX,self.curY))) 
        secondTerm = self.integral()
        logLike = firstTerm - secondTerm
        logLike = logLike + self.N + self.NlogDx + self.NlogDy
        mlogLike = logLike * (-1)
        #print("Parameters: %s -> %s (%s,%s)" %(parameters,mlogLike,firstTerm,secondTerm))
        return mlogLike
        
    def minimize(self):
        self.minuit.migrad()
        return self.minuit.get_fmin()['fval']
            
    def draw(self,tmin,tmax,xmin,xmax,ymin,ymax,Nev=100):
        self.setRange(tmin,tmax,xmin,xmax,ymin,ymax,0,self.t.shape[0])
        
        x = np.linspace(xmin, xmax, Nev)
        y = np.linspace(ymin, ymax, Nev)
        X, Y = np.meshgrid(x, y)
        zs = np.array([self.evaluate(0,x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
        Z = zs.reshape(X.shape)
        plt.contour(X,Y,Z)
       
class IsoModel(UnbinnedModel):
    
    def __init__(self,*args,**kwargs):
        super(IsoModel,self).__init__(*args,**kwargs)
        
        self.wrapper = lambda norm: self.mloglike(norm)
        
        #Assume uniform distribution in space and time
        init_value = self.N
        
        self.minuit = Minuit(self.wrapper, 
                             norm=init_value, limit_norm = (init_value/100.0,init_value*100.0),
                             error_norm = init_value/20.0,
                             fix_norm=False, print_level=-1,
                             errordef=0.5)
        self.minuit.set_strategy(0)
        
    def setParameters(self,norm):
        self.norm = norm
    
    def evaluate(self,t,x,y):
        return self.norm * np.ones(t.shape[0]) * 1.0/self.dOmega
    
    def integral(self):
        return self.norm
        
class GaussianIsoModel(UnbinnedModel):
    
    def __init__(self,*args,**kwargs):
        super(GaussianIsoModel,self).__init__(*args,**kwargs)
        
        self.wrapper = lambda x0,y0,sigma,norm,const: self.mloglike(x0,y0,sigma,norm,const)
        
        init_value = self.N/2.0
        
        self.minuit = Minuit(self.wrapper, 
                             x0=(self.xmin + self.dx/2.0), limit_x0 = (self.xmin,self.xmax),
                             error_x0 = 1.0,
                             y0=(self.ymin + self.dy/2.0), limit_y0 = (self.ymin,self.ymax),
                             error_y0 = 1.0,
                             sigma=90.0, fix_sigma = True,
                             norm=self.N/2.0, limit_norm = (1e-5,self.N*10), error_norm = self.N/100.0,
                             const = init_value, limit_const = (init_value/1000,init_value*1000),
                             error_const = init_value/20.0,
                             fix_const=False, print_level=-1,
                             errordef=0.5)
        self.minuit.set_strategy(0)

    
    def setParameters(self,x0,y0,sigma,norm,const):
        self.x0 = x0
        self.y0 = y0
        self.sigma = sigma
        self.sigmasq = pow(sigma,2)
        self.norm = norm
        self.const = const
        self._normal = multivariate_normal([self.x0,self.y0],np.identity(2)*self.sigmasq)
    
    def evaluate(self,t,x,y):
        #Equivalent to vstack, but faster
        reform = np.concatenate([[x,],[y,]],axis=0).T
        return self.norm/self.dt * self._normal.pdf(reform) + self.const * 1.0/self.dOmega
    
    def _integral(self,z0,a,b):
        #Integral of 1d gaussian between a and b
        sq2 = self.sigma * math.sqrt(2)
        one = -0.5 * erf( (z0 - a) / (sq2) )
        two = -0.5 * erf( (z0 - b) / (sq2) )
        return (two - one)
    
    def integral(self,):        
        return self._npredIso() + self._npredGauss()
    
    def _npredIso(self):
        return self.const
    
    def _npredGauss(self):
        xintegral = self._integral(self.x0, self.xmin, self.xmax)
        yintegral = self._integral(self.y0, self.ymin, self.ymax)
        return xintegral * yintegral * self.norm

"(1 / (1 + ( (x - xpos) / r0 )**2)**alpha) * (1 / (1 + ( (y - ypos) / r0 )**2)**alpha)"

class KingIsoModel(UnbinnedModel):
    
    def __init__(self,*args,**kwargs):
        super(KingIsoModel,self).__init__(*args,**kwargs)
        
        self.wrapper = lambda x0,y0,coreRadius,alpha,norm,const: self.mloglike(x0,y0,coreRadius,alpha,norm,const)
        
        init_value = self.N/2.0
        init_value = 1e-12
        
        self.minuit = Minuit(self.wrapper, 
                             x0=(self.xmin + self.dx/2.0), limit_x0 = (self.xmin,self.xmax),
                             error_x0 = 1.0, fix_x0 = False,
                             y0=(self.ymin + self.dy/2.0), limit_y0 = (self.ymin,self.ymax),
                             error_y0 = 1.0, fix_y0 = False,
                             coreRadius=85.0, fix_coreRadius = True, limit_coreRadius = (50,300),
                             alpha=1.5, fix_alpha=True, limit_alpha=(0,4),error_alpha = 0.05,
                             norm=self.N/2.0, limit_norm = (1e-5,self.N*1000.0), error_norm = self.N/100.0,
                             const = init_value, limit_const = (init_value/1000,init_value*1000),
                             error_const = init_value/20.0,
                             fix_const=True, print_level=-1,
                             errordef=0.5)
        self.minuit.set_strategy(0)

    
    def setParameters(self,x0,y0,coreRadius,alpha,norm,const):
        self.x0 = x0
        self.y0 = y0
        self.coreRadius = coreRadius
        self.alpha = alpha
        self.norm = norm
        self.const = const
        #This is the integral between -inf and inf of the king function,
        #in rdr:
        self.kingNorm = np.pi * pow(self.coreRadius,2.0) / (self.alpha - 1)
    
    def evaluate(self,t,x,y):        
        return self.norm/self.dt * self.evalKing(x,y) + self.const * 1.0/self.dOmega
    
    def evalKing(self,x,y):
        r = np.sqrt(np.power(x - self.x0,2.0) + np.power(y - self.y0,2.0))
        
        kingPDF = ( np.power(1 + np.power( r / self.coreRadius,2), (-1) * self.alpha) /
                  self.kingNorm )
        
        return kingPDF
        
    def integral(self,):        
        return self._npredIso() + self._npredKing()
    
    def _npredIso(self):
        return self.const
    
    #def _integral(self,rr):
    #    #King function integral on rdr
    #    return ( np.pi * pow(self.coreRadius,2.0) * 
    #             np.power(1 - np.power(rr/self.coreRadius,2),1 - self.alpha ) /
    #             (self.alpha - 1.0)  )

    
    def _npredKing(self):
        #Since there is no guarantee that the king function is all contained in the
        #region, I have to integrate numerically even though in principle the king
        #function is analytically integrable in rdr
        return self.norm * dblquad(self.evalKing,
                       self.xmin,self.xmax,
                       lambda x:self.ymin,lambda x:self.ymax)[0]

class SpaceClusteringLikelihood(object):
  def __init__(self,t,tmin,tmax,x,xmin,xmax,y,ymin,ymax,model=GaussianIsoModel):
    self.x                  = x
    self.y                  = y    
    self.t                  = t
    
    self.tmin               = tmin
    self.tmax               = tmax
    self.xmin               = xmin
    self.xmax               = xmax
    self.ymin               = ymin
    self.ymax               = ymax
    
    self.model              = model(t,tmin,tmax,x,xmin,xmax,y,ymin,ymax)
  pass
  
  def setPrior(self,ncp):
    self.ncp                = ncp
  
  def getPriors(self,N,p0):
    #return 4 - np.log(73.53 * p0 * np.power(np.arange(1,N+1),-0.478))
    return np.ones(N) * self.ncp
  
  def __call__(self,R,edges):
    
    Nevt = R+1
    
    allLikelihoods               = np.zeros(Nevt)
    TSs                          = np.zeros(Nevt)
    edg                          = []
    
    for i in range(Nevt):
      
      if( (Nevt) - i <= 0):
      
        thisNnLogL             = 1e6
        TS                     = 0
      else:
        self.model.setRange(edges[i],edges[R+1],
                            self.xmin,self.xmax,
                            self.ymin,self.ymax,
                            i,Nevt)
                
        thisLogL               = (-1) * self.model.minimize()
        TS                     = 2 * (thisLogL - getIsoLikelihoodValue(Nevt-i,edges[i],edges[R+1]))
        
      pass
      
      edg.append([edges[i],edges[R+1]])
      allLikelihoods[i]        = thisLogL
      TSs[i]                   = TS
      
    pass
    self.allLikelihoods          = allLikelihoods
    
    return allLikelihoods, TSs, edg
pass

