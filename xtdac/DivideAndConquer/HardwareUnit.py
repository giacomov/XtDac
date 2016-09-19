# Author: Giacomo Vianello (giacomov@stanford.edu)

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import logging
import subprocess

from xtwp4.DivideAndConquer import XMMWCS

log = logging.getLogger("HardwareUnit")


def hardwareUnitFactory(eventfile):
    '''Return an instance of the appropriate HarwareUnit class,
    based on the content of the event file.
    '''
    # Probe which hardware unit we are talking about
    # instance the corresponding class and return it

    with pyfits.open(eventfile) as f:
        instrume = f['EVENTS'].header.get("INSTRUME")
        if (instrume == "EPN"):
            return PNQuadrant(eventfile)
        elif (instrume == "EMOS1" or instrume == "EMOS2"):
            return MOSCCD(eventfile)
        elif instrume == "ACIS":
            return ACISCCD(eventfile)
        else:
            raise NotImplemented("Don't know how to handle instrument %s" % (instrume))
        pass


class HardwareUnit(object):
    '''
    A hardware unit is a part of a detector sharing one set of good time
    intervals, and roughly the same background. For example, a quadrant of the PN
    detector or one single CCD of one of the MOS are hardware units.
    '''

    def __init__(self, name, pixelScale,
                 minX, maxX,
                 minY, maxY):
        # Verify that all inputs are integers

        assert isinstance(minX, int)
        assert isinstance(maxX, int)

        assert isinstance(minY, int)
        assert isinstance(maxY, int)

        self.name = name

        self.minX = minX
        self.minY = minY
        self.maxX = maxX
        self.maxY = maxY

        self.pixelScale = float(pixelScale)  # arcsec

    def getName(self):
        return self.name

    def getMaxX(self):
        return self.maxX

    def getMinX(self):
        return self.minX

    def getMaxY(self):
        return self.maxY

    def getMinY(self):
        return self.minY

    def getPixelScale(self):
        return self.pixelScale


class ACISCCD(HardwareUnit):

    @staticmethod
    def get_ccd_name(id):
        # From http://cxc.harvard.edu/contrib/jcm/ncoords.ps
        ccd_ids = {0: 'I0',
                   1: 'I1',
                   2: 'I2',
                   3: 'I3',
                   4: 'S0',
                   5: 'S1',
                   6: 'S2',
                   7: 'S3',
                   8: 'S4',
                   9: 'S5'}

        return ccd_ids[id]

    def __init__(self, eventFile):
        # Open the event file
        with pyfits.open(eventFile) as f:
            # probe which CCD we are talking about

            ccd_id = f['EVENTS'].data.field("ccd_id")

            minccd, maxccd = ccd_id.min(), ccd_id.max()

            if minccd != maxccd:
                raise RuntimeError("Provided event file contains events from more than one ACIS CCD.")

            if minccd > 9:
                raise RuntimeError("The provided event file is not a ACIS event file")

            X = f['EVENTS'].data.field("X")
            Y = f['EVENTS'].data.field("Y")

        name = "ACIS (CCD %s)" % self.get_ccd_name(minccd)
        self.ccd_id = minccd

        # ACIS pixels are 0.492 arcsec

        super(ACISCCD, self).__init__(name, 0.492,
                                      int(X.min()), int(X.max()),
                                      int(Y.min()), int(Y.max()))

    def get_psf(self, xpos, ypos, eventfile, outfile):
        """
        Produce an image containing the PSF for the given position

        :return: the filename containing the image
        """

        # First run the MARX simulation to get the image

        # Get R.A., Dec. of the source

        wcs = XMMWCS.XMMWCS(eventfile)

        ra, dec = wcs.xy2sky([[xpos, ypos]])[0]

        # Figure out the name of the detector (ACIS-I or ACIS-S)
        # CCDs 0 to 3 are ACIS-I, the others (up to 9) are ACIS-S
        if self.ccd_id <= 3:

            detector_type = 'ACIS-I'

        else:

            detector_type = 'ACIS-S'

        # Get the pointing
        ra_pnt = wcs.ra_nom
        dec_pnt = wcs.dec_nom
        rot_pnt = wcs.rotation

        # Now get the position of the SIM (http://cxc.harvard.edu/ciao/threads/marx_sim/)

        header = pyfits.getheader(eventfile, "EVENTS")

        obs_sim_x = header['SIM_X']
        obs_sim_y = header['SIM_Y']
        obs_sim_z = header['SIM_Z']

        if detector_type == 'ACIS-I':

            marx_sim_x, marx_sim_y, marx_sim_z = -0.7823481983384, 0, -233.5924630914

        else:

            marx_sim_x, marx_sim_y, marx_sim_z = -0.68426746699586, 0, -190.1325231040

        delta_sim_x = obs_sim_x - marx_sim_x
        delta_sim_y = obs_sim_y - marx_sim_y
        delta_sim_z = obs_sim_z - marx_sim_z

        # Get the DY and DZ from the eventfile
        try:

            dx = 0
            dy = header['DY_AVG']
            dz = header['DZ_AVG']

        except KeyError:

            raise RuntimeError("DY and DZ keywords not found in %s. You have to run r4_header_update on it.")

        # Finally we can compute the offsets for MARX

        DetOffsetX = dx + delta_sim_x
        DetOffsetY = dy + delta_sim_y
        DetOffsetZ = dz + delta_sim_z

        # Run the simulation detecting 100000 events (using a negative NumRays means that the simulation
        # will continue until 100000 events are detected)

        cmd_line = "marx ExposureTime=0.0 NumRays=-10000 GratingType=NONE SourceRA=%s SourceDEC=%s MinEnergy=1.5 " \
                   "MaxEnergy=1.5 SourceFlux=1 SourceType=POINT OutputDir=__point DetectorType=%s " \
                   "RA_Nom=%s Dec_Nom=%s Roll_Nom=%s DetOffsetX=%s DetOffsetY=%s DetOffsetZ=%s Verbose=no" \
                   % (ra, dec, detector_type, ra_pnt, dec_pnt, rot_pnt, DetOffsetX, DetOffsetY, DetOffsetZ)

        log.debug(cmd_line)

        subprocess.check_call(cmd_line, shell=True)

        # Now generate the FITS file

        cmd_line = 'marx2fits __point __sim.fits'
        log.debug(cmd_line)

        _ = subprocess.check_output(cmd_line, shell=True)

        # Finally generate the FITS image
        xmin, xmax = xpos - 500, xpos + 500
        ymin, ymax = ypos - 500, ypos + 500

        cmd_line = "f2dhisto __sim.fits %s 1 1 X Y '%s,%s' '%s,%s' clobber=yes" % (outfile, xmin, xmax, ymin, ymax)
        log.debug(cmd_line)
        _ = subprocess.check_output(cmd_line, shell=True)

        return outfile


class XMMCCD(object):

    @staticmethod
    def get_psf(xpos, ypos, eventfile, outfile):
        """
        Produce an image containing the PSF for the given position

        :return: the filename containing the image
        """

        pars = {}
        pars['output'] = outfile
        pars['region'] = '"(X,Y) IN circle(%s,%s)"' % (xpos, ypos)
        pars['image'] = eventfile
        pars['xsize'] = 500
        pars['ysize'] = 500
        pars['energy'] = '"200 600 1500 6000 10000"'
        pars['level'] = 'ELLBETA'
        pars['-V'] = 0

        cmdLine = 'psfgen %s' % (" ".join(['%s=%s' % (k, v) for k, v in pars.iteritems()]))

        log.debug(cmdLine)
        subprocess.check_call(cmdLine, shell=True)

        return outfile


class PNQuadrant(HardwareUnit, XMMCCD):
    def __init__(self, eventFile):
        quadrant = None

        # Open the event file

        with pyfits.open(eventFile, memmap=False) as f:

            # probe which quadrant we are talking about
            ccdnr = f['EVENTS'].data.field("CCDNR")
            (minccd, maxccd) = (ccdnr.min(), ccdnr.max())
            if (minccd >= 1 and maxccd <= 3):
                # Quadrant 1
                quadrant = 1
                (mindetx, maxdetx) = (-18283, -2260)
                (mindety, maxdety) = (-1090, 15325)
            elif (minccd >= 4 and maxccd <= 6):
                # Quadrant 2
                quadrant = 2
                (mindetx, maxdetx) = (-2150, 13880)
                (mindety, maxdety) = (-1090, 15325)
            elif (minccd >= 7 and maxccd <= 9):
                # Quadrant 3
                quadrant = 3
                (mindetx, maxdetx) = (-2150, 13880)
                (mindety, maxdety) = (-17527, -1110)
            elif (minccd >= 10 and maxccd <= 12):
                # Quadrant 4
                quadrant = 4
                (mindetx, maxdetx) = (-18283, -2260)
                (mindety, maxdety) = (-17527, -1110)
            else:
                raise RuntimeError(
                    "The provided event file %s contains events from more than one PN quadrant." % (eventFile))

            X = f['EVENTS'].data.field("X")
            Y = f['EVENTS'].data.field("Y")


        name = "EPIC PN (quadrant %i)" % quadrant

        super(PNQuadrant, self).__init__(name, 0.05,
                                         int(X.min()), int(X.max()),
                                         int(Y.min()), int(Y.max())
                                         )


class MOSCCD(HardwareUnit, XMMCCD):
    def __init__(self, eventFile):
        ccdNumber = None

        # Open the event file
        with pyfits.open(eventFile) as f:

            # Get the CCD number
            ccdnr = f['EVENTS'].data.field("CCDNR")

            # Verify that there is only one CCD number in the event file
            # otherwise throw an exception

            uniqueCCDnumbers = set(ccdnr)

            if (len(uniqueCCDnumbers) > 1):

                raise RuntimeError("The provided event file contains events from more than one MOS CCD!")

            else:

                ccdNumber = ccdnr[0]

            pass

            # Read the event coordinates

            X = f['EVENTS'].data.field("X")
            Y = f['EVENTS'].data.field("Y")

        pass

        name = "EPIC MOS (CCD %i)" % ccdNumber

        super(MOSCCD, self).__init__(name, 0.05,
                                     int(X.min()), int(X.max()),
                                     int(Y.min()), int(Y.max())
                                     )
