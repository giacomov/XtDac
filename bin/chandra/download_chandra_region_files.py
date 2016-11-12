#!/usr/bin/env python

import subprocess
import numpy as np
import os
import sys
import warnings
import multiprocessing
import time


def worker(obsid):

    cmd_line = 'obsid_search_csc obsid=%s outfile=_.tsv download=all ' \
               'filetype=reg mode=h clobber=yes verbose=0' % obsid

    if os.path.exists(str(obsid)):

        return obsid

    for j in range(10):

        try:

            subprocess.check_call(cmd_line, shell=True)

        except subprocess.CalledProcessError:

            time.sleep(5.0)

            continue

        else:

            break

    if not os.path.exists(str(obsid)):


        return None

    else:

        return obsid

if __name__=="__main__":

    # Read observation ids

    obsids = np.recfromtxt("all_obs.txt", names=True)

    # Create output directory (if doesn't exist)

    if not os.path.exists('region_files'):

        os.makedirs('region_files')

    orig_dir = os.getcwd()

    os.chdir('region_files')

    pool = multiprocessing.Pool(8)

    good = []
    bad = []

    for i, res in enumerate(pool.imap(worker, obsids['obsid'])):

        if res is not None:

            good.append(res)

        else:

            bad.append(obsids['obsid'][i])

        sys.stderr.write("\r%s out of %s (%s bad)" % (i + 1, len(obsids['obsid']), len(bad)))

    os.chdir(orig_dir)

    with open("good_obsids.txt","w+") as f:

        for g in good:

            f.write("%s\n" % g)


