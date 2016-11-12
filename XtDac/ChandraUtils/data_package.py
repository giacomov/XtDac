import fnmatch
import os
import shutil

import yaml

from XtDac.ChandraUtils.logging_system import get_logger
from XtDac.ChandraUtils.sanitize_filename import sanitize_filename
from XtDac.ChandraUtils.work_within_directory import work_within_directory

logger = get_logger("DataPackage")

_index_file = "index.yml"


def _check_directory(directory):
    sanitized_directory = sanitize_filename(directory)

    assert os.path.exists(sanitized_directory), "Directory %s does not exists" % sanitized_directory

    assert os.path.isdir(sanitized_directory), "The file %s is not a directory" % sanitized_directory

    return sanitized_directory


class File(object):
    def __init__(self, filename, description):
        self._filename = sanitize_filename(filename)
        self._description = description

        assert os.path.exists(self._filename), "Something went wrong when creating File instance. " \
                                               "File %s does not exists!" % self._filename

    @property
    def filename(self):
        return self._filename

    @property
    def description(self):
        return self._description

    def _check_consistency(self):
        assert os.path.exists(self._filename), "The file %s was moved independently, bad idea!" % self._filename

    def move_to(self, new_directory):
        """
        Move the file to a new location (a new directory)

        :param new_directory:
        :return: new path
        """

        self._check_consistency()

        new_directory = _check_directory(new_directory)

        shutil.move(self._filename, new_directory)

        new_path = os.path.join(new_directory, os.path.basename(self._filename))

        assert os.path.exists(new_path), "Could not move %s to %s" % (self._filename, new_path)

        self._filename = new_path

        return self

    def copy_to(self, new_directory):
        self._check_consistency()

        new_directory = _check_directory(new_directory)

        shutil.copy(self._filename, new_directory)

        new_path = os.path.join(new_directory, os.path.basename(self._filename))

        assert os.path.exists(new_path), "Could not copy %s to %s" % (self._filename, new_path)

        return File(new_path, self._description)


class DataPackage(object):
    def __init__(self, directory, create=False):

        self._directory = sanitize_filename(directory)

        if os.path.exists(self._directory) and os.path.isdir(self._directory):

            logger.debug("Accessing data in %s" % self._directory)

            with work_within_directory(self._directory):

                # Access the index file
                assert os.path.exists(_index_file), "Cannot find index file in %s" % self._directory

                self._load_status()

                self._check_consistency()

        else:

            if create:

                # Create directory

                os.makedirs(self._directory)

                # Create an empty index file

                with work_within_directory(self._directory):

                    # By default the package is read-write

                    self._status = {'read_only': False, 'index': {}}

                    self._save_status()

                logger.info("Datapackage in %s has been created" % self._directory)

            else:

                raise IOError("Directory %s does not exist or is not a directory" % self._directory)

    @property
    def location(self):

        return self._directory

    def has(self, tag):
        """
        Returns whether the package contains the file corresponding to the provided tag, or not

        :param tag:
        :return: True or False
        """

        self._load_status()

        return tag in self._status['index']

    def _set_readonly(self, read_only):

        read_only = bool(read_only)

        self._status['read_only'] = read_only

        self._save_status()

    def _get_readonly(self):

        return self._status['read_only']

    read_only = property(_get_readonly, _set_readonly, doc="Set or get the read-only status of the package")

    @property
    def _status_file(self):

        return os.path.abspath(os.path.join(self._directory, _index_file))

    def _get_abs_path(self, tag):

        self._load_status()

        return os.path.abspath(os.path.join(self._directory, self._status['index'][tag]['path']))

    def _save_status(self):

        # Save the dictionary to the dictionary file

        with open(self._status_file, "w+") as f:
            yaml.dump(self._status, f)

    def _load_status(self):

        # Read the dictionary

        with open(self._status_file, "r") as f:
            self._status = yaml.load(f)

    def _check_consistency(self):

        self._load_status()

        # Check that all files described in the dictionary exist

        with work_within_directory(self._directory):

            for tag in self._status['index'].keys():

                path = self._status['index'][tag]['path']

                if not os.path.exists(path):
                    abspath = os.path.abspath(path)

                    raise IOError("File %s is contained in the index, but does not exists in %s" % (path, abspath))

    def clear(self):
        """
        Remove all files from the data package (careful!)

        :return: None
        """

        if self.read_only:
            raise RuntimeError("Trying to modifying a read-only package")

        self._check_consistency()

        for tag in self._status['index'].keys():
            # NOTE: the order here is important, because _get_abs_path reload the status!

            path = self._get_abs_path(tag)

            self._status['index'].pop(tag)

            os.remove(path)

            self._save_status()

    def store(self, tag, filename, description, force=False, move=False):
        """
        Store (move) a file in the package

        :param filename:
        :param description:
        :return:
        """

        self._load_status()

        if self.read_only:
            raise RuntimeError("Trying to modifying a read-only package")

        if tag in self._status['index'] and not force:
            raise RuntimeError("Cannot store file with tag %s, because the tag is already present in the package. "
                               "Use .update()." % tag)

        # Create the instance of a File

        orig_file = File(filename, description)

        # Move the file inside the package

        if move:

            new_file = orig_file.move_to(self._directory)

        else:

            new_file = orig_file.copy_to(self._directory)

        # Register it in the dictionary (using a relative path)

        relative_path = os.path.relpath(new_file.filename, self._directory)

        self._status['index'][tag] = {'path': relative_path, 'description': orig_file.description}

        # Save to the index file

        self._save_status()

    def update(self, tag, filename):
        """
        Update a file which is already in the package

        :param tag:
        :param filename:
        :return:
        """

        self._load_status()

        if self.read_only:
            raise RuntimeError("Trying to modifying a read-only package")

        assert tag in self._status['index'], "Cannot update file with tag %s, it does not exist in the package" % tag

        # Assure that the filename is the same

        name1 = os.path.basename(filename)
        name2 = os.path.basename(self._status['index'][tag]['path'])

        if name1 != name2:
            raise RuntimeError("You cannot update the file %s with tag %s with the file %s "
                               "which has a different name" % (name1, tag, name2))

        # Move old file to a temporary location

        temp_backup = os.path.join(self._directory, name1 + '.bak')

        path = self._get_abs_path(tag)

        shutil.move(path, temp_backup)

        # Store new one

        try:

            self.store(tag, filename, self._status['index'][tag]['description'], force=True)

        except:

            # Move back the temp file
            shutil.move(temp_backup, path)

            logger.error("Could not update file with tag %s, could not store the new file. "
                         "The old file has been restored." % tag)

            raise

        else:

            # If we are here the store has worked out fine, remove the temp file
            os.remove(temp_backup)

    def get(self, tag, dest_dir=None):
        """
        Retrieve a file by tag from the data package

        :param tag:
        :param dest_dir: if None, use current workdir, otherwise use the one provided, as destination dir. for the file
        :return: a File instance
        """

        self._load_status()

        assert tag in self._status['index'], "Tag %s does not exists in data package: \n%s" % (tag, self)

        item = self._status['index'][tag]

        abs_path = self._get_abs_path(tag)

        this_file = File(abs_path, item['description'])

        if dest_dir is not None:

            dest = sanitize_filename(dest_dir)

        else:

            dest = os.getcwd()

        out_file = this_file.copy_to(dest)

        return out_file

    def copy_to(self, new_directory):
        """
        Copy the entire data package to another directory

        :param new_directory: destination path. The package will be moved, with its name, inside this directory,
        which must already exist
        :return: the instance of the new package
        """

        self._save_status()

        directory = _check_directory(new_directory)

        package_name = os.path.split(self._directory)[-1]

        destination = os.path.join(directory, package_name)

        shutil.copytree(self._directory, destination)

        return DataPackage(destination)

    def __repr__(self):

        self._load_status()

        repr = ""
        repr += "Data package in %s" % self._directory

        if self._status['read_only']:

            repr += " (read only)\n"

        else:

            repr += " (read/write)\n"

        for tag in self._status['index']:
            repr += "* %s: %s\n" % (tag, self._status['index'][tag]['path'])

        return repr

    def find_all(self, pattern):
        """
        Returns all tags matching the patter

        :param pattern: a pattern like "ccd_*" (unix-style wildcards)
        :return: list of tags matching the pattern
        """

        tags = fnmatch.filter(self._status['index'].keys(), pattern)

        return tags
