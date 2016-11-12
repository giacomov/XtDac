#!/usr/bin/env python

""""
Take event file and create multiple new event files separated by CCD

command from CIAO: dmcopy filtered_event.fits[EVENTS][ccd_id=N] out.fits clobber=yes

Make sure CIAO is running before running this script
"""

import argparse
import os
import subprocess

import astropy.io.fits as pyfits
from XtDac.ChandraUtils.run_command import CommandRunner
from XtDac.ChandraUtils.logging_system import get_logger

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create new event files separated by CCD')

    parser.add_argument('--evtfile', help="Event file name", type=str, required=True)

    args = parser.parse_args()

    logger = get_logger("xtc_separate_CCDs.py")
    runner = CommandRunner(logger)

    print "Separating by CCD..."

    for ccd_id in range(10):

        ccd_file = "ccd_%s_%s" % (ccd_id, os.path.basename(args.evtfile))

        cmd_line = "dmcopy %s[EVENTS][ccd_id=%s] %s clobber=yes" % (args.evtfile, ccd_id, ccd_file)

        runner.run(cmd_line)

        # check if certain CCD files are empty and then delete them if so

        f = pyfits.open("%s" % (ccd_file))
        ccd_data = f[1].data

        if len(ccd_data) == 0:
            os.remove(ccd_file)

        f.close()
