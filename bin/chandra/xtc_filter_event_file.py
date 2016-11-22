#!/usr/bin/env python

"""
Take evt3 file and use region files to subtract off sources that are already known - image will have lots of holes

Make sure CIAO is running before running this script

below = code used by Giacomo to create filtered image
ftcopy 'acisf00635_000N001_evt3.fits[EVENTS][regfilter("my_source.reg")]' test.fits

code that works with CIAO
dmcopy "acisf00635_000N001_evt3.fits[exclude sky=region(acisf00635_000N001_r0101_reg3.fits)]" filter_test.fits opt=all
"""

import argparse
import glob
import os
import shutil
import re
import sys
import warnings

import astropy.io.fits as pyfits
import numpy as np

from XtDac.ChandraUtils import find_files
from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils import sanitize_filename
from XtDac.ChandraUtils import setup_ftools
from XtDac.ChandraUtils.data_package import DataPackage
from XtDac.ChandraUtils.run_command import CommandRunner


def is_variable(tsv_file, name_of_the_source):
    # Read data from tsv file
    tsv_data = np.array(np.recfromtxt(tsv_file, delimiter='\t', skip_header=11, names=True), ndmin=1)

    # Make sure data is vectorized
    if len(tsv_data.shape) == 0:
        tsv_data = np.array([tsv_data])

    tsv_names = tsv_data['name'].tolist()

    # Find index of source in question
    idx = tsv_names.index(name_of_the_source)

    # See if associated var_flag is True (some are strings with spaces, some are bools)
    if str(tsv_data['var_flag'][idx]) == "TRUE" or str(tsv_data['var_flag'][idx]) == " TRUE":
        return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Filter known sources out of event file')

    parser.add_argument('--in_package', help="Data package directory for the input data", type=str, required=True)

    parser.add_argument('--region_dir', help="Directory containing the regions file for this obsid",
                        type=str, required=True)

    parser.add_argument('--out_package', help="Data package directory for the output data", type=str, required=True)

    parser.add_argument("--debug", help="Debug mode (yes or no)", type=bool, required=False, default=False)

    parser.add_argument("--emin", help="Minimum energy (eV)", type=int, required=True)

    parser.add_argument("--emax", help="Maximum energy (eV)", type=int, required=True)

    parser.add_argument("--adj_factor",
                        help="If region files need to be adjusted, what factor to increase axes of ellipses by",
                        type=float, required=True)

    parser.add_argument("--randomize_time", help='Whether to randomize the arrival time within the time frames',
                        required=True, dest='randomize_time', action='store_true')

    parser.set_defaults(randomize_time=True)

    # assumption = all level 3 region files and event file are already downloaded into same directory

    # Get logger for this command

    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    # Instance the command runner

    runner = CommandRunner(logger)

    args = parser.parse_args()

    # Setup the FTOOLS so they can be run non-interactively

    setup_ftools.setup_ftools_non_interactive()

    # creates text file with name of all level 3 region files for given Obs ID

    region_dir = sanitize_filename.sanitize_filename(args.region_dir)

    obsid = os.path.split(region_dir)[-1]

    # region_dir is specific to one obsid. Get general region repository where db is located

    db_dir = os.path.split(region_dir)[0]

    # Get the region files from this observation

    region_files = find_files.find_files(region_dir, "*reg3.fits.gz")

    # Open the data package
    data_package = DataPackage(args.in_package)

    # Get the pointing from the event file

    evtfile = data_package.get('evt3').filename
    fovfile = data_package.get('fov3').filename

    with pyfits.open(evtfile) as f:

        ra_pnt = f['EVENTS'].header.get("RA_PNT")
        dec_pnt = f['EVENTS'].header.get("DEC_PNT")

        frame_time = f['EVENTS'].header['TIMEDEL']

        logger.info("Found a frame time of %s" % (frame_time))

    n_reg = len(region_files)

    logger.info("Found %s sources in the CSC for this obsid (%s)" % (n_reg, obsid))

    # Loop over the region files and prepare the corresponding region for filtering

    temp_files = []

    # This will be set to True if there is at least one source which might have produced streaks of out-of-time events
    might_have_streaks = False

    temp_filter = "__temp_filter_evt.fits"

    shutil.copy2(evtfile, temp_filter)

    for region_id, region_file in enumerate(region_files):

        if region_id % 50 == 0 or region_id == len(region_files) - 1:
            sys.stderr.write("\rProcessing region %s out of %s ..." % (region_id + 1, len(region_files)))

        temp_outfile = '__temp_filter_evt_post.fits'

        cmd_line = 'dmcopy \"%s[exclude sky=region(%s)][opt update=no]\" ' \
                   '%s opt=all clobber=yes' % (temp_filter, region_file, temp_outfile)

        runner.run(cmd_line, debug=True)

        os.remove(temp_filter)
        os.rename(temp_outfile, temp_filter)

        # Copy the region file here and fix the column format, otherwise ftcopy later on will fail
        # because some files have vector X,Y columns, while some others don't

        temp_file = "__%i_%s_reg.fits" % (region_id, obsid)

        shutil.copy2(region_file, "%s.gz" % temp_file)

        # Unzip the file because fcollen cannot operate on zipped files
        if os.path.exists(temp_file):

            os.remove(temp_file)

        runner.run("gunzip %s.gz" % temp_file, debug=True)

        # Fix the columns

        cmd_line = "fcollen '%s' X 1" % temp_file
        runner.run(cmd_line, debug=True)

        cmd_line = "fcollen '%s' Y 1" % temp_file
        runner.run(cmd_line, debug=True)

        temp_files.append(temp_file)

        ##############$$$######################
        # Check if it might have caused streaks
        #################$$$###################

        # Avoid checking again if we know that there are already streaks (the presence of one streaking source
        # will cause the run of the de-streaking tool whether or not there are other streaking sources)

        if not might_have_streaks:

            # Check if the source is on a streak. If it is, it means that this observation has streaks (simple ah?)
            with pyfits.open(region_file) as f:

                if bool(f['SRCREG'].header.get('ONSTREAK')) == True:

                    might_have_streaks = True


    sys.stderr.write("Done\n")

    # Write all region files in a text file, which will be used to create a merged file
    # with all regions
    regions_list_file = "__all_regions_list.txt"

    with open(regions_list_file, "w+") as f:

        for region in temp_files:
            f.write("%s\n" % region)

    # Merge all the region files

    logger.info("Merging region files...")

    all_regions_file = '%s_all_regions.fits' % (obsid)

    cmd_line = 'ftmerge @%s %s clobber=yes columns=-' % (regions_list_file, all_regions_file)

    runner.run(cmd_line)

    # Now fix the COMPONENT column (each region must be a different component)

    fits_file = pyfits.open(all_regions_file, mode='update', memmap=False)

    fits_file['SRCREG'].data.COMPONENT[:] = range(fits_file['SRCREG'].data.shape[0])

    fits_file.close()

    # Add the region file to the output package

    output_package = DataPackage(args.out_package, create=True)

    output_package.store('all_regions', all_regions_file, "FITS file containing the regions relative to all sources "
                                                          "for this OBSID")

    # Finally filter the file

    logger.info("Filtering event file...")

    ###########################
    # Filter by energy
    ###########################

    outfile = '%s_filtered_evt3.fits' % obsid

    cmd_line = 'dmcopy %s[energy=%d:%d] %s opt=all clobber=yes' % (temp_filter, args.emin, args.emax, outfile)

    runner.run(cmd_line)

    ###########################
    # Remove readout streaks
    ###########################

    # If we might have streaks, run the tool which generates the region for the streaks

    if might_have_streaks:

        logger.warn("We might have readout streak. Cleaning them out...")

        streak_region_fileroot = "streakreg"

        streak_region_ds9 = "%s.reg" % streak_region_fileroot

        # Run the tool which produces a region file with the streak
        cmd_line = 'acis_streak_map infile=%s fovfile=%s bkgroot=__bkg ' \
                   'regfile=%s clobber=yes msigma=5' % (outfile, fovfile, streak_region_ds9)

        runner.run(cmd_line)

        # Check if acis_streak_map has found anything
        with open(streak_region_ds9) as f:

            lines = f.readlines()

        if len(filter(lambda x: x.find("Polygon") >= 0, lines)) == 0:

            # No streak found
            logger.warn("No streak found by acis_streak_map")

        else:

            # Transform in FITS format

            streak_region_file = "%s.fits" % (streak_region_fileroot)

            cmd_line = "dmmakereg 'region(%s.reg)' %s kernel=fits clobber=yes" % (streak_region_fileroot,
                                                                                  streak_region_file)

            runner.run(cmd_line)

            # Now change the name of the extension to conform to the other region files

            with pyfits.open(streak_region_file) as f:

                header = f['REGION'].header

                header.set("EXTNAME", "SRCREG")
                header.set("HDUNAME", "SRCREG")

                f.writeto(streak_region_file, clobber=True)

            temp_out = '__temp_out.fits'

            cmd_line = 'dmcopy \"%s[exclude sky=region(%s)]\" ' \
                       '%s opt=all clobber=yes' % (outfile, streak_region_file, temp_out)

            runner.run(cmd_line)

            os.remove(outfile)

            os.rename(temp_out, outfile)

            # Store the streak regions in the output package
            output_package.store("streak_regions_ds9", streak_region_ds9,
                                 "Regions containing streaks from bright point sources (out-of-time events)")

    ##############################
    # Randomize time             #
    ##############################

    with pyfits.open(outfile, mode='update') as f:

        time = f['EVENTS'].data.time

        logger.info("Randomizing arrival times within time frameof %s s..." % frame_time)

        deltas = np.random.uniform(-frame_time/2.0, frame_time/2.0, time.shape[0])

        time += deltas

        f["EVENTS"].data.time[:] = time

    # Now sort the file
    cmd_line = "fsort %s[EVENTS] TIME heap" % outfile

    runner.run(cmd_line)

    # Store output in the output package
    output_package.store("filtered_evt3", outfile, "Event file (Level 3) with all point sources in the CSC, "
                                                   "as well as out-of-time streaks (if any), removed")

    # remove files
    files_to_remove = glob.glob("__*")

    if not args.debug:

        for file in files_to_remove:

            try:

                os.remove(file)

            except:

                logger.warn("Could not remove %s" % file)



    else:

        logger.warn("\n\nWARNING: did not remove temporary files because we are in debug mode")
