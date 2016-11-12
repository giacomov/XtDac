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
import re
import sys
import warnings

import astropy.io.fits as pyfits
import numpy as np

from XtDac.ChandraUtils import find_files
from XtDac.ChandraUtils import logging_system
from XtDac.ChandraUtils import query_region_db
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


def cross_match(region_files_db, region_files_obsid):
    # Get the names of all the sources contained respectively in the OBSID catalog and in the DB catalog

    names_obsid = map(lambda x: re.match('(.+)/(CXOJ.+)/.+', x).groups()[1], region_files_obsid)

    names_db = map(lambda x: re.match('(.+)/(CXOJ.+)/.+', x).groups()[1], region_files_db)

    # Loop over all sources in the OBSID catalog, and remove them from teh DB catalog. At the end,
    # the union between the OBSID catalog and the remaining of the DB catalog will give the
    # complete catalog

    for source in names_obsid:

        try:

            index = names_db.index(source)

            if args.debug:
                print names_db[index]

            names_db.pop(index)

        except ValueError:

            warnings.warn("Source %s is in the OBSID catalog but not in the DB catalog" % source)

        else:

            if args.debug:
                print region_files_db[index]

            region_files_db.pop(index)

    cleaned_regions = list(region_files_obsid)

    cleaned_regions.extend(region_files_db)

    return cleaned_regions


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

    region_files_obsid = find_files.find_files(region_dir, "*reg3.fits.gz")

    # Open the data package
    data_package = DataPackage(args.in_package)

    # Get the pointing from the event file

    evtfile = data_package.get('evt3').filename
    fovfile = data_package.get('fov3').filename

    with pyfits.open(evtfile) as f:

        ra_pnt = f['EVENTS'].header.get("RA_PNT")
        dec_pnt = f['EVENTS'].header.get("DEC_PNT")

    # Query a region of 30 arcmin, which should always cover the whole Chandra field of view,
    # to get the regions from the database

    print "Querying region database..."

    region_files_db = query_region_db.query_region_db(ra_pnt, dec_pnt, 30.0, db_dir)

    # Now cross match the regions we got from the DB with the regions we got from this obsid
    # We try to use the information relative to this obsid as much as possible, but if there is no
    # info on a given source in this obsid we take it from the db

    print "Cross-matching region files..."

    region_files = cross_match(region_files_db, region_files_obsid)

    n_reg = len(region_files)

    # Loop over the region files and prepare the corresponding region for filtering

    temp_files = []

    # This will be set to True if there is at least one source which might have produced streaks of out-of-time events
    might_have_streaks = False

    for region_id, region_file in enumerate(region_files):

        if region_id % 50 == 0 or region_id == len(region_files) - 1:
            sys.stderr.write("\rProcessing region %s out of %s ..." % (region_id + 1, len(region_files)))

        temp_file = '__%d_reg_revised.fits' % (region_id)

        # Remove the file if existing
        if os.path.exists(temp_file):
            os.remove(temp_file)

        cmd_line = '''ftcopy '%s[SRCREG][SHAPE=="Ellipse"]' %s clobber=yes''' % (region_file, temp_file)

        runner.run(cmd_line, debug=True)

        # Fix the column format, if needed

        cmd_line = "fcollen '%s' X 1" % temp_file
        runner.run(cmd_line, debug=True)

        cmd_line = "fcollen '%s' Y 1" % temp_file
        runner.run(cmd_line, debug=True)

        # Adjust the size of the ellipse, if this source is variable

        # Get the name of the source
        # Filenames are absolute paths like
        # /home/giacomov/science/chandra_transients/catalog/region_files/635/CXOJ162659.0-243557/[filename]
        # so the second-last element is the name of the source

        source_name = os.path.basename(os.path.split(region_file)[-2])
        source_name = source_name[0:3] + " " + source_name[3::]

        ##########################
        # Crossmatch with tsv file
        ##########################

        if args.adj_factor > 1:

            try:

                if is_variable(data_package.get('tsv').filename, source_name) == True:
                    # open the file with "mode='update'"

                    with pyfits.open(temp_file, mode='update') as reg:
                        reg['SRCREG'].data.R = args.adj_factor * reg['SRCREG'].data.R

                        # data, h = fitsio.read(temp_file, ext='SRCREG', header=True)
                        #
                        # data['R']  = args.adj_factor * data['R']
                        #
                        # fitsio.write(temp_file, data, extname='SRCREG', header=h, clobber=True)

                        # adjust the size of both axis by a factor (another argument)

            except ValueError:

                pass

        temp_files.append(temp_file)

        ##############$$$######################
        # Check if it might have caused streaks
        #################$$$###################

        # Avoid checking again if we know that there are already streaks (the presence of one streaking source
        # will cause the run of the de-streaking tool whether or not there are other streaking sources)

        if not might_have_streaks:

            # Open the region file and get the countrate for this source

            with pyfits.open(temp_file) as f:

                # Total number of counts in the region
                roi_cnts = f['SRCREG'].header.get("ROI_CNTS")

                # Exposure
                exposure = f['SRCREG'].header.get("EXPOSURE")

                # Get the average rate
                rate = roi_cnts / float(exposure)

                # The readout time of the whole array is 3.2 s

                target_rate = 1.0 / 3.2

                if rate >= target_rate / 10.0:
                    # Source might generate out-of-time (readout streak) events

                    might_have_streaks = True

    sys.stderr.write("Done\n")

    # Write all the temp files in a text file, which will be used by dmmerge
    regions_list_file = "__all_regions_list.txt"

    with open(regions_list_file, "w+") as f:

        for temp_file in temp_files:
            f.write("%s\n" % temp_file)

    # Merge all the region files

    print("Merging region files...")

    all_regions_file = '%s_all_regions.fits' % (obsid)

    cmd_line = 'ftmerge @%s %s clobber=yes columns=-' % (regions_list_file, all_regions_file)

    runner.run(cmd_line)

    # Now fix the COMPONENT column (each region must be a different component, otherwise
    # dmcopy will crash)

    fits_file = pyfits.open(all_regions_file, mode='update', memmap=False)

    fits_file['SRCREG'].data.COMPONENT[:] = range(fits_file['SRCREG'].data.shape[0])

    fits_file.close()

    # Add the region file to the output package

    output_package = DataPackage(args.out_package)

    output_package.store('all_regions', all_regions_file, "FITS file containing the regions relative to all sources "
                                                          "for this OBSID")

    # Finally filter the file

    print("Filtering event file...")

    ###########################
    # Filter by energy
    ###########################

    temp_filter = '__temp__%s' % obsid

    cmd_line = 'dmcopy %s[energy=%d:%d] %s opt=all clobber=yes' % (evtfile, args.emin, args.emax, temp_filter)

    runner.run(cmd_line)

    ###########################
    # Filter by regions
    ###########################

    outfile = '%s_filtered_evt3.fits' % obsid

    cmd_line = 'dmcopy \"%s[exclude sky=region(%s)]\" ' \
               '%s opt=all clobber=yes' % (temp_filter, all_regions_file, outfile)

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
                   'regfile=%s clobber=yes msigma=3' % (outfile, fovfile, streak_region_ds9)

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

        frame_time = f['EVENTS'].header['TIMEDEL']

        logger.info("Found a frame time of %s. Randomizing arrival times within time frame..." % (frame_time))

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

        print("\n\nWARNING: did not remove temporary files because we are in debug mode")
