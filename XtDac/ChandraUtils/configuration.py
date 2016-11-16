import yaml
import os
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename


class ReadOnlyContainer(object):

    def __init__(self, dictionary):

        self._dict = dict(dictionary)

    def __getitem__(self, item):

        return self._dict[item]


def get_configuration(filename):

    filename = sanitize_filename(filename)

    assert os.path.exists(filename), "Configuration file %s does not exist!" % filename

    try:

        with open(filename, "r") as f:

            configuration_dict = yaml.safe_load(f)

    except:

        raise IOError("Couldn't read configuration file %s. File is not readable, or wrong format." % (filename))

    configuration = ReadOnlyContainer(configuration_dict)

    return configuration