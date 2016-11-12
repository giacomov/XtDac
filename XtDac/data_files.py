import pkg_resources
import os


def get_data_file_path(data_file):
    """
    Returns the absolute path to the required data files.

    :param data_file: relative path to the data file, relative to the XtDac/data path.
    So to get the path to data/test.dat you need to use data_file="test.dat"
    :return: absolute path of the data file
    """

    try:

        file_path = pkg_resources.resource_filename("XtDac", 'data/%s' % data_file)

    except KeyError:

        raise IOError("Could not read or find data file %s. Try reinstalling XtDac. If this does not fix your "
                      "problem, open an issue on github." % (data_file))

    else:

        return os.path.abspath(file_path)
