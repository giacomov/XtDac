import iminuit
import math
import matplotlib.pyplot as plt
import numexpr
import numpy
import scipy.interpolate

try:

    import nlopt

except:

    # We can live without
    pass

try:

    import scipy.optimize

except:

    # We can live without
    pass

import cartesian

import logging

import collections

from xtwp4.DivideAndConquer import XMMWCS

import ExposureMap

try:

    from astropy.io import fits as pyfits

except:

    import pyfits

log = logging.getLogger("main")


class Parameter(object):
    def __init__(self, name, value, error, minval, maxval, free=True):
        self.name = str(name)
        self.value = float(value)
        self.error = float(error)
        self.minval = minval
        self.maxval = maxval

        self.free = bool(free)

    def isFree(self):
        return self.free

    def fix(self):
        self.free = False

    def free(self):
        self.free = True

        # def __eq__(self, value):

        #    self.value = float( value )


class Model(object):
    def __init__(self, name, *parameters):

        self.name = str(name)

        self.parameters = collections.OrderedDict()

        for p in parameters:
            self.parameters[p.name] = p

    def getName(self):

        return self.name

    def getParameters(self):

        return self.parameters

    def setParameters(self, values):

        assert len(values) == len(self.parameters.keys()), \
            "Trying to set parameters with incorrect number of values"

        for k, v in zip(self.parameters.keys(), values):
            self[k] = v

    def __setitem__(self, key, value):

        self._verifyKey(key)

        self.parameters[key].value = value

    def _verifyKey(self, key):

        assert key in self.parameters.keys(), \
            "Trying to access parameter %s which is not part of the current model" % key

    def __getitem__(self, key):

        self._verifyKey(key)

        return self.parameters[key]


class Isotropic(Model):
    def __init__(self, name, value):
        p = Parameter("%s_amp" % (name), value, value / 10.0, 1e-6, None)

        super(Isotropic, self).__init__(name, p)

    def __call__(self, x, y):
        out = numpy.zeros((x.shape[0], y.shape[0]), dtype=float)

        out.fill(self["%s_amp" % (self.name)].value)

        return out.T


class PointSource(Model):

    def __init__(self, name, xpos, ypos, eventfile, hwu, psffile=None):

        # NOTE: we assume that xpos and ypos cannot change

        self.xpos = float(xpos)
        self.ypos = float(ypos)

        if psffile is None:

            # Generate the PSF image at this position

            outfile = "%s_PSF.fits" % name

            hwu.get_psf(xpos, ypos, eventfile, outfile)

        else:

            outfile = psffile

        # Now read the FITS file
        with pyfits.open(outfile) as f:

            psfimg = f[0].data

            xsize = psfimg.shape[1]
            ysize = psfimg.shape[0]

        # Generate the interpolator

        # Transform the coordinates from image to X,Y

        pixScale = hwu.getPixelScale()

        minx = self.xpos - (xsize / 2.0 / pixScale)
        maxx = self.xpos + (xsize / 2.0 / pixScale)

        miny = self.ypos - (ysize / 2.0 / pixScale)
        maxy = self.ypos + (ysize / 2.0 / pixScale)

        xx = numpy.linspace(minx, maxx, xsize)
        yy = numpy.linspace(miny, maxy, ysize)

        self.psf = scipy.interpolate.interp2d(xx, yy, psfimg)

        # Now initiate the mother class
        amp = Parameter("%s_norm" % name, 1.0, 0.1, 0, None)

        self.outfile = outfile

        super(PointSource, self).__init__(name, amp)

    def __call__(self, x, y):

        out = self["%s_norm" % self.name].value * self.psf(x, y)

        return out


class GlobalModel(object):
    def __init__(self, name):

        self.name = str(name)

        self.sources = collections.OrderedDict()

    def __add__(self, src):

        self.addSource(src)

        return self

    def addSource(self, src):

        assert isinstance(src, Model), "The new source must be a Model instance!"

        assert src.name not in self.sources.keys(), "Source %s is already present!" % src.name

        self.sources[src.getName()] = src

        self._updateParameters()

    def removeSource(self, srcName):

        assert srcName in self.sources.keys(), "Source %s is not part of this global model." % srcName

        src = self.sources.pop(srcName)

        self._updateParameters()

        return src

    def _updateParameters(self):

        self.parameters = collections.OrderedDict()

        for k, src in self.sources.iteritems():
            p = src.getParameters()

            self.parameters.update(p)

    def __getitem__(self, paramName):

        assert paramName in self.parameters.keys(), \
            "Parameter %s does not belong to this model" % paramName

        return self.parameters[paramName]

    def getParameterNames(self):

        return self.parameters.keys()

    def setParameters(self, parameters):

        assert len(parameters) == len(self.parameters.keys()), \
            "wrong number of parameters in setParameters()"

        for key, value in zip(self.parameters.keys(), parameters):
            self.parameters[key].value = value

    def __call__(self, x, y):

        out = numpy.sum([src(x, y) for src in self.sources.values()], axis=0)

        return out


class Region(object):
    def __init__(self, xmin, xmax, ymin, ymax, binsize, regfilter, expomap, eventfile):

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        assert self.ymax > self.ymin
        assert self.xmax > self.xmin

        # At least two bins!
        assert binsize < (self.xmax - self.xmin) / 2.0
        assert binsize < (self.ymax - self.ymin) / 2.0

        self.binsize = binsize

        self.xbins = numpy.arange(self.xmin, self.xmax, self.binsize)
        self.ybins = numpy.arange(self.ymin, self.ymax, self.binsize)

        # Prepare mask

        self.xcs = (self.xbins[:-1] + self.xbins[1:]) / 2.0
        self.ycs = (self.ybins[:-1] + self.ybins[1:]) / 2.0

        points = cartesian.cartesian([self.xcs, self.ycs])

        # plt.imshow(exposures)

        self.mask = numpy.zeros((self.ycs.shape[0], self.xcs.shape[0]), dtype=bool)

        for i in range(self.mask.shape[0]):

            for j in range(self.mask.shape[1]):
                self.mask[i, j] = regfilter.inside1(self.xcs[j], self.ycs[i])

        expo = ExposureMap.Exposure(expomap)

        wcs = XMMWCS.XMMWCS(eventfile)

        skypoints = wcs.xy2sky(points)
        ras = skypoints[:, 0]
        decs = skypoints[:, 1]
        exposures = expo.getExposureAtSky(ras, decs).reshape(self.mask.T.shape).T

        idx = exposures <= (exposures.max() / 2.0)

        self.mask[idx] = 0

        self.empty = ~self.mask


class Likelihood(object):
    def __init__(self, x, y, region):

        self.region = region

        self.img, _, _ = numpy.histogram2d(x, y, [region.xbins, region.ybins])
        self.img = self.img.T

        self.like = numpy.zeros_like(self.img)

        assert self.img.shape == self.region.mask.shape, \
            "Img shape: %s, mask shape: %s" % (self.img.shape, self.region.mask.shape)

    def getTotalCounts(self):

        return numpy.sum(self.img[self.region.mask])

    def getMaximumPosition(self):

        pixy, pixx = numpy.unravel_index(self.img.argmax(), self.img.shape)

        x = self.region.xcs[pixx]
        y = self.region.ycs[pixy]

        return x, y

    def setModel(self, lmodel):

        # Model must be LikelihoodModel
        assert isinstance(lmodel, GlobalModel), \
            "Model must be an instance of GlobalModel"

        self.model = lmodel

    def mloglike(self, *values):

        self.model.setParameters(values)

        self.m = self.model(self.region.xcs, self.region.ycs)

        assert self.m.shape == self.img.shape

        # Binned likelihood

        # o * log(m) - m - log(o!)

        im = self.img[self.region.mask]
        mm = self.m[self.region.mask]

        # lf = gammaln( self.img[ self.region.mask ] + 1 )

        ll = numexpr.evaluate("im * log( mm ) - mm")

        return numpy.sum(ll) * (-1)

    def minimize(self, minuit=True, verbose=1):

        if minuit:

            paramNames = self.model.getParameterNames()

            # Build Minuit parameters specification
            self.mpars = {}
            self.mpars['forced_parameters'] = paramNames
            self.mpars['errordef'] = 0.5
            self.mpars['print_level'] = verbose

            for k in self.mpars['forced_parameters']:
                self.mpars[k] = self.model[k].value
                self.mpars["error_%s" % k] = self.model[k].error
                self.mpars["limit_%s" % k] = [self.model[k].minval, self.model[k].maxval]
                self.mpars["fix_%s" % k] = (not self.model[k].free)

            minuit = iminuit.Minuit(self.mloglike, **self.mpars)

            minuit.tol = 1000
            minuit.set_strategy(0)

            minuit.migrad()

            # Update the parameters with the results of the fit
            for k in paramNames:
                self.model[k].value = minuit.values[k]
                self.model[k].error = minuit.errors[k]

            return minuit.fval

        elif 1 == 0:

            newpars = list(self.model.parameters)

            def wrap(args):

                newpars[0] = args[0]
                newpars[1] = args[1]

                logl = self.mloglike(*newpars)

                # print("%s -> %s" %(args, logl))

                return logl

            res = scipy.optimize.minimize(wrap, [0.1, 1.0], method='L-BFGS-B',
                                          bounds=[[0, 10], [0, 10]])

            # res = scipy.optimize.minimize_scafor you to do it on my behalf ?lar( wrap, bounds=(0, 10), method='bounded')

            return res.fun

        else:

            newpars = list(self.model.parameters)

            def wrapper(args, grad):

                if grad.size > 0:
                    print("This won't ever happen, since BOBYQA does not use derivatives")

                newpars[0] = args[0]
                newpars[1] = args[1]

                logl = self.mloglike(*newpars)

                return logl

            bob = nlopt.opt(nlopt.LN_NELDERMEAD, 2)

            bob.set_ftol_abs(0.3)

            # Stop if the value of all the parameter change by less than 1%
            # bob.set_xtol_rel(1e-2)
            bob.set_initial_step([0.1, 0.1])

            bob.set_lower_bounds([0, 0])
            bob.set_upper_bounds([10, 10])

            bob.set_min_objective(wrapper)

            res = bob.optimize([0.1, 1.0])
            opt_val = bob.last_optimum_value()
            result = bob.last_optimize_result()

            return opt_val

    def plot(self):

        vmin = self.img.min()
        vmax = self.img.max()

        fig = plt.figure(figsize=(12, 5))
        plt.subplot(1, 3, 1)
        self.img[self.region.empty] = numpy.nan
        plt.imshow(self.img, interpolation='nearest', vmin=vmin, vmax=vmax, origin='lower')
        plt.title("Data")

        plt.subplot(1, 3, 2)
        self.m[self.region.empty] = numpy.nan
        plt.imshow(self.m, interpolation='nearest', vmin=vmin, vmax=vmax, origin='lower')
        plt.title("Model")

        plt.subplot(1, 3, 3)
        res = (self.img - self.m) / numpy.sqrt(self.m)
        res[self.region.empty] = numpy.nan
        plt.imshow(res, interpolation='nearest', vmin=-3, vmax=3, origin='lower')
        plt.colorbar()
        plt.title("Residuals")

        return fig


class LikelihoodScanner(object):
    def __init__(self, x, xmin, xmax, y, ymin, ymax, filter=None):

        self.x = numpy.asarray(x)
        self.y = numpy.asarray(y)

        assert self.x.shape[0] == self.y.shape[0]

        if filter is None:

            self.l = Likelihood(x, xmin, xmax, y, ymin, ymax)

        else:

            self.l = Likelihood(x, xmin, xmax, y, ymin, ymax, filter=filter)

        # x1,y1 = self.l.getPixelPos(22848.00753125, 26592.56254865)
        # x2,y2 = self.l.getPixelPos(23252.32863361, 26795.65978228)



        self.model = LikelihoodModel(self.l._x.min(), self.l._x.max(),
                                     self.l._y.min(), self.l._y.max())

        # self.model.addKnownSource( x1, y1 )
        # self.model.addKnownSource( x2, y2 )

        log.info("\nLikelihood scanner initiated")
        log.debug("Total counts: %s" % (self.x.shape[0]))
        log.debug("Counts surviving filter: %s" % (self.l.getTotalCounts()))
        log.debug("Median in image: %s\n" % (numpy.median(self.l.img[~self.l.empty].flatten())))

        # plt.imshow(self.l.empty)
        # plt.imshow(self.l.img)

        self.knownSources = 0

    def addKnownSource(self, X, Y, *args):

        x, y = self.l.getPixelPos(X, Y)
        self.model.addKnownSource(x, y, *args)
        self.knownSources += 1

        log.debug("Known source position: (%s,%s)" % (x, y))

    def nullModelFit(self):

        self.model.removeTestSource()
        self.l.setModel(self.model)
        like0 = self.l.maximize(minuit=True)

        return like0

    def getBackgroundLevel(self):

        like0 = self.nullModelFit()

        value = self.model.model.amplitude_0.value / pow(float(self.l.binsize), 2)
        error = self.l.minuitErrors['amplitude_0'] / pow(float(self.l.binsize), 2)

        return value, error, like0

    def findExcesses(self, nx=5, ny=5):

        xpoints = numpy.linspace(self.l._x.min(), self.l._x.max(), nx)
        ypoints = numpy.linspace(self.l._y.min(), self.l._y.max(), ny)

        like0 = self.nullModelFit()

        # Fix known sources
        for i in range(self.knownSources):
            thisAmpName = 'amplitude_%s' % (self.knownSources + 1)
            self.model.mpars[thisAmpName] = self.model.model.__getattr__(thisAmpName).value
            self.model.mpars['fix_amplitude_%s' % (self.knownSources + 1)] = True

        self.model.addTestSource()

        bestPos = [0, 0]
        tsmap = numpy.zeros([nx, ny])

        for i, x in enumerate(xpoints):

            for j, y in enumerate(ypoints):

                d = math.sqrt(pow(x - self.l.xc, 2) + pow(y - self.l.yc, 2))

                if d > self.l.xc:
                    continue

                self.model.setTestSourcePosition(x, y)

                self.l.setModel(self.model)
                like1 = self.l.maximize()
                xc, yc = (self.model.model.x_mean_1.value, self.model.model.y_mean_1.value)

                thisTS = 2 * (like0 - like1)

                # print("bkg = %s" %(self.model.model.amplitude_0.value))
                # print("pos = (%s, %s) -> TS = %s" %(xc, yc, thisTS))

                # print("(X,Y) = (%s,%s) -> TS = %s" %(x,y,thisTS))

                tsmap[i, j] = thisTS

        log.info("Max TS = %s" % (tsmap.max()))

        i, j = numpy.unravel_index(tsmap.argmax(), tsmap.shape)
        # print i,j

        xsize = (xpoints[1] - xpoints[0])
        ysize = (ypoints[1] - ypoints[0])

        xxc, yyc = xpoints[i] + xsize / 2.0, ypoints[j] + ysize / 2.0

        self.centroids = numpy.array([self.l.getPosOfPixel(xxc, yyc)])

        if 1 == 1:
            plt.figure(figsize=(12, 5))
            plt.subplot(1, 2, 1)
            plt.imshow(tsmap, interpolation='nearest')
            plt.title("TS map")
            plt.subplot(1, 2, 2)
            plt.imshow(self.l.img, interpolation='nearest')
            plt.plot([yyc], [xxc], 'o', color='red')
            plt.title("Image")

        return self.centroids

    def plot(self, xbins, ybins):

        fig, sub = plt.subplots(1, 1)

        _ = sub.hist2d(self.x, self.y, [self.l.xbins, self.l.ybins])

        if len(self.centroids) > 0:
            cc = numpy.array(self.centroids)

            sub.plot(cc[:, 0], cc[:, 1], 'o', markersize=15, color='red')

        return fig
