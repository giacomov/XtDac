import os

import numpy as np

from XtDac.ChandraUtils.angular_distance import angular_distance


def query_region_db(ra_center, dec_center, radius, region_dir):
    """
    Returns a list of files relative to regions which are within the provided cone

    :param ra_center: R.A. of the center of the cone
    :param dec_center: Dec. of the center of the cone
    :param radius: Radius of the cone (in arcmin)
    :param region_dir: Path of the directory containing the database file and the region files
    :return: list of region files
    """

    database_file = os.path.join(region_dir, 'region_database.txt')

    data = np.recfromtxt(database_file, names=True)

    # Compute the angular distance between all regions and the center of the cone

    distances = angular_distance(float(ra_center), float(dec_center), data["RA"], data["DEC"], unit='arcmin')

    # Select all regions within the cone

    idx = (distances <= float(radius))

    # Return the corresponding region files

    data_list = []

    for i in data['REGION_FILE'][idx]:
        data_list.append(os.path.join(region_dir, i))

    return data_list
