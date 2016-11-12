import astropy.units as u
from astropy.coordinates import SkyCoord


def angular_distance(ra1, dec1, ra2, dec2, unit='degree'):
    """
    Compute angular distance between pairs of coordinates

    :param ra1:
    :param dec1:
    :param ra2:
    :param dec2:
    :param unit: the unit of the output (default: degree)
    :return: angular distance
    """
    point_1 = SkyCoord(ra=ra1 * u.degree, dec=dec1 * u.degree, frame='icrs')

    point_2 = SkyCoord(ra=ra2 * u.degree, dec=dec2 * u.degree, frame='icrs')

    angular_distance = point_1.separation(point_2).to(unit)

    return angular_distance.value
