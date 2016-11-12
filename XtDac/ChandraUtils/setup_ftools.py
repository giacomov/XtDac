import os


def setup_ftools_non_interactive():
    """
    Setup the environment so the FTOOLS can be run in the batch farm without a terminal attached to them.

    :return: none
    """

    os.environ['HEADASNOQUERY'] = ''
    os.environ['HEADASPROMPT'] = ''
