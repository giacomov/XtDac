import fnmatch
import os


def find_files(directory, pattern):
    matches = []

    for root, dirnames, filenames in os.walk(directory):

        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.abspath(os.path.join(root, filename)))

    return matches
