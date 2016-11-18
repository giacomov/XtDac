import subprocess

from XtDac.ChandraUtils.setup_ftools import setup_ftools_non_interactive

class CommandRunner(object):
    def __init__(self, logger):

        self._logger = logger

    def run(self, cmd_line, debug=False):

        # Make sure the FTOOLS are non-interactive, otherwise some commands will fail in the computer farm
        setup_ftools_non_interactive()

        if debug:

            self._logger.debug(cmd_line)

        else:

            self._logger.info(cmd_line)

        subprocess.check_call(cmd_line, shell=True)
