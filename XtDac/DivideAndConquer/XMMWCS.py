# Author: Giacomo Vianello (giacomov@stanford.edu)

from astropy import wcs as pywcs

try:
    import astropy.io.fits as pyfits
except:
    # If this fail there is no way out
    import pyfits
pass
import numpy
from Box import Box


class XMMWCS(object):

    def __init__(self, eventFile, X=None, Y=None):

        header = pyfits.getheader(eventFile, 'EVENTS')

        self.wcs = pywcs.WCS(naxis=2)

        if 'REFXCRPX' in header:

            # XMM

            self.wcs.wcs.crpix = [header['REFXCRPX'], header['REFYCRPX']]
            self.wcs.wcs.cdelt = [header['REFXCDLT'], header['REFYCDLT']]
            self.wcs.wcs.crval = [header['REFXCRVL'], header['REFYCRVL']]
            self.wcs.wcs.ctype = [header['REFXCTYP'], header['REFYCTYP']]

            # Position angle is the angle between the R.A. axis and the horizontal axis
            # of the CCD
            self.rotation = header['PA_PNT']



        else:

            # Chandra
            try:

                self.wcs.wcs.crpix = [header['TCRPX11'], header['TCRPX12']]
                self.wcs.wcs.cdelt = [header['TCDLT11'], header['TCDLT12']]
                self.wcs.wcs.crval = [header['TCRVL11'], header['TCRVL12']]
                self.wcs.wcs.ctype = [header['TCTYP11'], header['TCTYP12']]

            except KeyError:

                # Try harder, by looking first where the 'x' and 'y' columns are
                # (in a simulated file they could be in a different position)
                x_pos = int(filter(lambda x: str(x[1]).replace(" ", "") == "x",
                                   header.items())[0][0].replace("TTYPE", ""))
                y_pos = int(filter(lambda x: str(x[1]).replace(" ", "") == "y",
                                   header.items())[0][0].replace("TTYPE", ""))

                self.wcs.wcs.crpix = [header['TCRPX%i' % x_pos], header['TCRPX%i' % y_pos]]
                self.wcs.wcs.cdelt = [header['TCDLT%i' % x_pos], header['TCDLT%i' % y_pos]]
                self.wcs.wcs.crval = [header['TCRVL%i' % x_pos], header['TCRVL%i' % y_pos]]
                self.wcs.wcs.ctype = [header['TCTYP%i' % x_pos], header['TCTYP%i' % y_pos]]

            self.rotation = header['ROLL_PNT']

        self.ra_nom = header['RA_PNT']
        self.dec_nom = header['DEC_PNT']

        if X is None or Y is None:

            with pyfits.open(eventFile) as f:

                X = f['EVENTS'].data.field("X")
                Y = f['EVENTS'].data.field("Y")

        self.xmin = X.min()
        self.xmax = X.max()
        self.ymin = Y.min()
        self.ymax = Y.max()

    pass

    def xy2sky(self, points):
        '''
      Convert from X-Y coordinates to R.A., Dec.

      Usage:
        ra, dec = XMMWCS.xy2sky(points)
      where:
        points                    : list of X,Y pairs. For example: points = [[1000,2000],[3000,5000]]
                                    represents two points (1000,2000) and (3000,5000)
      Returns:
        a numpy array with the same number of R.A.,Dec. pairs as the input points array
        '''
        pixcrd = numpy.array(points, numpy.float_)
        sky = self.wcs.all_pix2world(pixcrd, 1)
        return sky

    pass

    def sky2xy(self, points):
        '''
      Convert from R.A.,Dec coordinates to X-Y coordinates

      Usage:
        x, y = XMMWCS.sky2xy(points)
      where:
        points                    : list of Ra,Dec pairs. For example: points = [[10.72,-9.75],[13.25,179.5]]
                                    represents two points (10.72,-9.75) and (13.25,179.5)
      Returns:
        a numpy array with the same number of X,Y pairs as the input points array
        '''
        pixcrd = numpy.array(points, numpy.float_)
        sky = self.wcs.wcs_world2pix(pixcrd, 1)
        return sky

    pass

    def makeGrid(self, side=80, outfile=None):
        '''
        Make a grid in R.A., Dec. made of boxes with the given side (in arcsec). The
        grid will be parallel to the axes of the CCD.
        '''
        # Get coordinates of the center
        ra_c, dec_c = self.wcs.wcs.crval

        # Now build the grid aligned with the R.A.,Dec. axis
        # Start from the center, and move of 'side' step
        #    boxes                     = []
        #    PNside                    = 30.0 #arcmin
        #    Nside                     = int(numpy.ceil(PNside*60.0/float(side)))
        #    #Nside                     = 1
        #
        #    #sky                       = self.xy2sky([[self.xmin,self.ymin]])
        #    #ra_corner,dec_corner      = sky[0]
        #
        #    ra_corner                 = ra_c - (PNside/60.0/2.0)
        #    dec_corner                = dec_c - (PNside/60.0/2.0)
        #
        #    for i in range(Nside):
        #      for j in range(Nside):
        #        boxes.append(Box(ra_corner+i*(side/3600.0),dec_corner+j*(side/3600.0),side,side,0,'fk5'))
        #      pass
        #    pass

        side_pix = side / 0.05

        xs = numpy.arange(self.xmin, self.xmax, side_pix / 2.0)
        ys = numpy.arange(self.ymin, self.ymax, side_pix / 2.0)
        grid_ = numpy.meshgrid(xs, ys)
        centers_pix = zip(*(x.flat for x in grid_))

        regionsDs9 = []
        boxes = []

        for i, (x, y) in enumerate(centers_pix):
            box = Box(x, y, side / 0.05, side / 0.05, 0, 'physical')
            boxes.append(box)
            regionsDs9.append(box.getDs9Region()[-1])

        pass

        if (outfile != None):
            with open(outfile, "w+") as f:
                f.write("physical\n")
                f.write("\n".join(regionsDs9))

        return boxes

    pass


pass
