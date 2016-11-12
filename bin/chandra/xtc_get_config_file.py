#!/usr/bin/env python

"""
Copy the configuration file from the repository to the current directory, for editing.
"""

import argparse
import os
import sys

from XtDac.ChandraUtils import logging_system
from XtDac.data_files import get_data_file_path

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Get a copy of the template for the configuration file')

    parser.add_argument("-o", "--output", help="Name for the new configuration file", type=str, required=True)

    # Get the logger
    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    args = parser.parse_args()

    # Load configuration file from the code repository
    copy_of_the_file = []

    with open(get_data_file_path('sample_configuration.yml')) as f:

        template_lines = f.readlines()

    # Dump the sample dict in the output file
    # NOTE: we do not use a copy because the file in the code repository might be (likely will be) read-only, while
    # we want to create a file that the user can actually edit
    with open(args.output, 'w+') as f:

        f.writelines(template_lines)
