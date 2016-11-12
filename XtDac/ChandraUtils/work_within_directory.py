import contextlib
import os


@contextlib.contextmanager
def work_within_directory(directory):
    original_directory = os.getcwd()

    os.chdir(directory)

    try:

        yield

    except:

        raise

    finally:

        os.chdir(original_directory)
