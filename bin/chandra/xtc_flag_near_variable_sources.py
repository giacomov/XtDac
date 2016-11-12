#!/usr/bin/env python

"""
Take output from Bayesian block algorithm and cross check matches with previously flagged variable X-ray sources

NOTE: There may be more than one obsid listed for the same source
"""

import argparse
import os
import sys

import astropy.units as u
import numpy as np

from XtDac.ChandraUtils import chandra_psf
from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils import offaxis_angle
from XtDac.ChandraUtils.chandra_catalog import ChandraSourceCatalog
from XtDac.ChandraUtils.run_command import CommandRunner

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Take output from Bayesian block algorithm and cross match with "
                                                 "previously flagged variable sources")
    parser.add_argument("--bbfile", help="Path of input text file (checked for hot pixels)", required=True)
    parser.add_argument("--outfile", help="Path of out file", required=True)
    parser.add_argument("--eventfile", help="Event file (needed to gather the pointing of Chandra)", required=True,
                        type=str)

    # Get logger for this command

    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Instance the command runner

    runner = CommandRunner(logger)

    args = parser.parse_args()

    # Instance the catalog
    csc = ChandraSourceCatalog()

    # get directory path and file name from input file arguments

    bb_file_path = os.path.abspath(os.path.expandvars(os.path.expanduser(args.bbfile)))

    # read BB data into array
    bb_data = np.array(np.recfromtxt(bb_file_path, names=True), ndmin=1)

    # number of rows of data
    bb_n = len(bb_data)

    # Create a PSF instance, we will use it to compute the size of the PSF at the position
    # of the source

    psf = chandra_psf.ChandraPSF()

    with open(args.outfile, "w") as f:

        # Pre-existing column names
        existing_column_names = " ".join(bb_data.dtype.names)

        f.write("# %s Closest_Variable_Source Separation(arcsec) Var_msid Theta PSF_size(arcsec) PSFfrac\n"
                % existing_column_names)

        for i in xrange(bb_n):

            # Compute the radius of the PSF at the off-axis angle of this source

            hotpix_flag = bb_data['Hot_Pixel_Flag'][i]
            ra = bb_data['RA'][i]
            dec = bb_data['Dec'][i]

            theta = offaxis_angle.get_offaxis_angle(ra, dec, args.eventfile)  # arcmin

            psf_size = psf.get_psf_size(theta, percent_level=0.95)

            if hotpix_flag == True:

                temp_list = []

                for j in xrange(len(bb_data.dtype.names)):
                    temp_list.append(str(bb_data[i][j]))

                # Fill the columns "Closest_Variable_Source","Separation","Var_msid", "Theta", "PSF","PSFfrac"
                # with appropriate info

                temp_list.append("None")
                temp_list.append(str(-1))
                temp_list.append(str(0))
                temp_list.append(str(theta))
                temp_list.append(str(psf_size))
                temp_list.append(str(1))

                line = " ".join(temp_list)

                f.write("%s\n" % line)

            else:

                # search_csc has a max search radius of 60, so put upper bound on input

                radius = 5.0

                variable_sources = csc.find_variable_sources(ra, dec, radius, unit='arcmin', column='var_flag')

                # Find source with minimum separation

                if variable_sources.shape[0] == 0:

                    # No variable source in the catalog. Let's add a "None" for each entry in the new columns

                    temp_list = []

                    for j in xrange(len(bb_data.dtype.names)):
                        temp_list.append(str(bb_data[i][j]))

                    # Fill the columns "Closest_Variable_Source","Separation","Var_msid", "Theta", "PSF","PSFfrac"
                    # with appropriate info

                    temp_list.append("None")
                    temp_list.append(str(-1))
                    temp_list.append(str(0))
                    temp_list.append(str(theta))
                    temp_list.append(str(psf_size))
                    temp_list.append(str(1))

                    line = " ".join(temp_list)

                    f.write("%s\n" % line)

                else:

                    closest_variable_source = csc.find_closest_variable_source(ra, dec)

                    # Get the name/separation/msid of the closest variable source

                    src_name = closest_variable_source['name']
                    src_sepn = (closest_variable_source['distance'] * u.arcmin).to(u.arcsec).value
                    src_msid = closest_variable_source['msid']

                    # Replace any space in the name with an underscore
                    src_name = src_name.replace(" ", "_")

                    psf_frac = psf.get_psf_fraction(theta, src_sepn)

                    temp_list = []

                    for j in xrange(len(bb_data.dtype.names)):
                        temp_list.append(str(bb_data[i][j]))

                    # Fill the columns "Closest_Variable_Source","Separation","Var_msid", "Theta,"PSF","PSFfrac"
                    # with appropriate info

                    temp_list.append(src_name)
                    temp_list.append(str(src_sepn))
                    temp_list.append(str(src_msid))
                    temp_list.append(str(theta))
                    temp_list.append(str(psf_size))
                    temp_list.append(str(psf_frac))

                    line = " ".join(temp_list)

                    f.write("%s\n" % line)
