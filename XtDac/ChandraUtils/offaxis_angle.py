# Chandra. Automatically compute the step size, as the average
# size of the PSF in the detector

import astropy.io.fits as pyfits
import astropy.units as u
from astropy.coordinates import SkyCoord


def get_offaxis_angle(ra, dec, event_file):
    # Get the aim point

    with pyfits.open(event_file) as f:
        ra_pointing = f['EVENTS'].header.get("RA_PNT")
        dec_pointing = f['EVENTS'].header.get("DEC_PNT")
        system = f['EVENTS'].header.get("RADECSYS")

    if system is None:
        system = 'ICRS'

    # Compute the corresponding off-axis angle theta

    pointing = SkyCoord(ra=ra_pointing * u.degree, dec=dec_pointing * u.degree, frame=system.lower())

    c1 = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame=system.lower())

    this_theta = c1.separation(pointing)

    # Store it in arcmin

    theta = this_theta.to(u.arcmin).value

    return theta
