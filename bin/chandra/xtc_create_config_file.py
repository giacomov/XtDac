#!/usr/bin/env python

"""
Copy the configuration file from the repository to the current directory, for editing.
"""

import argparse
import os
import sys
import yaml
import textwrap

from XtDac.ChandraUtils import logging_system
from XtDac.data_files import get_data_file_path
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create a configuration file')

    parser.add_argument("-o", "--output", help="Name for the new configuration file", type=str, required=True)

    # Now get the parameters from the config file and add them as command line options
    conf_file_template = get_data_file_path('sample_configuration.yml')

    with open(conf_file_template, "r") as f:

        tokens = []

        current_comment = ''
        current_name = ''
        current_value = ''

        for line in f.xreadlines():

            if line[0]=='\n':

                # Empty line, skip
                continue

            elif line[0]=='#':

                # This is a comment
                current_comment += line.replace("#"," ").replace("%", " percent")

            else:

                # Value
                current_name = line.split(":")[0].lstrip().rstrip().replace(" ","_")
                current_value = line.split(":")[1].lstrip().replace("\n","")

                tokens.append((current_comment.lstrip().rstrip().replace("\n"," ").replace("  "," "),
                               current_name, current_value))

                current_comment = ''
                current_name = ''
                current_value = ''

    for (comment, name, value) in tokens:

        if name.find("repository") > 0:

            parser.add_argument("--%s" % name, help=comment, required=True)

        else:

            parser.add_argument("--%s" % name, help=comment, default=value)

    # Get the logger
    logger = logging_system.get_logger(os.path.basename(sys.argv[0]))

    args = parser.parse_args()

    # Now make sure that all 3 repositories are different from each other
    # NOTE: this is a set, hence its entries are kept unique
    directories = {sanitize_filename(args.data_repository),
                   sanitize_filename(args.region_repository),
                   sanitize_filename(args.output_repository)}

    assert len(directories) == 3, "The data, region and output repositories must point at different directories"

    # Load configuration file from the code repository
    with open(conf_file_template) as f:

        template = yaml.safe_load(f)

    with open(args.output, "w+") as f:

        for (comment, name, value) in tokens:

            # Write comment

            comment_lines = textwrap.wrap(comment, 80)

            for line in comment_lines:

                f.write("# %s\n" % line)

            name_with_spaces = name.replace("_", " ")

            assert name_with_spaces in template

            # Write key : value pair

            if name_with_spaces.find("repository") > 0 or name_with_spaces == 'work directory':

                # Take absolute path
                abs_path = sanitize_filename(args.__getattribute__(name))

                # Create if needed
                if not os.path.exists(abs_path):

                    logger.info("Directory %s does not exist, creating it" % abs_path)
                    os.makedirs(abs_path)

                f.write("%s : %s\n\n" % (name_with_spaces, abs_path))

            else:

                f.write("%s : %s\n\n" % (name_with_spaces, args.__getattribute__(name)))