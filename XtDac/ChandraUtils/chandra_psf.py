import caldb4
import psf


class ChandraPSF(object):
    def __init__(self):
        cdb = caldb4.Caldb(telescope="CHANDRA", product="REEF")
        reef = cdb.search[0]
        extno = cdb.extno()

        # Replace the trailing '[..]' block number specifier
        reef = reef.split('[')[0] + "[{}]".format(extno + 1)

        self._pdata = psf.psfInit(reef)

    def get_psf_size(self, angle_in_arcmin, energy_in_kev=1.5, percent_level=0.9, phi=0.0):
        # Get the size of the psf at the minimum and maximum theta
        # for an energy of 1.5 keV

        psf_size = psf.psfSize(self._pdata, energy_in_kev, angle_in_arcmin, phi, percent_level)

        return psf_size

    # returns PSF in arcsec

    def get_psf_fraction(self, angle_in_arcmin, distance_in_arcsec, energy_in_kev=1.5, phi=0.0):
        return psf.psfFrac(self._pdata, energy_in_kev, angle_in_arcmin, phi, distance_in_arcsec)
