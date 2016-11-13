import contextlib
import os
import shutil
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename

@contextlib.contextmanager
def work_within_directory(directory, create=False, remove=False):

    directory = sanitize_filename(directory)

    original_directory = os.getcwd()

    if not os.path.exists(directory):

        if create:

            try:

                os.makedirs(directory)

            except:

                raise IOError("Cannot create directory %s" % directory)
        else:

            raise IOError("Directory %s does not exist" % directory)

    os.chdir(directory)

    try:

        yield

    except:

        raise

    finally:

        os.chdir(original_directory)

        if remove:

            shutil.rmtree(directory)