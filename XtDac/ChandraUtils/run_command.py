import subprocess


class CommandRunner(object):
    def __init__(self, logger):

        self._logger = logger

    def run(self, cmd_line, debug=False):

        if debug:

            self._logger.debug(cmd_line)

        else:

            self._logger.info(cmd_line)

        subprocess.check_call(cmd_line, shell=True)
