#!/usr/bin/env python

import glob
import os
from setuptools import setup

# Get the version number
execfile('XtDac/version.py')

# Normal packages

packages = ['XtDac',
            'XtDac/BayesianBlocks',
            'XtDac/DivideAndConquer',
            'XtDac/FixedBinSearch',
            'XtDac/ChandraUtils'
            ]

# XMM scripts
xmm_scripts = glob.glob(os.path.join('bin','xmm','*.py'))
# Chandra scripts
chandra_scripts = glob.glob(os.path.join('bin','chandra','*.py'))
# General executables
main_executables = glob.glob(os.path.join('bin','*.py'))

setup(
    name="XtDac",

    packages=packages,

    # data_files=[('astromodels/data/functions', glob.glob('astromodels/data/functions/*.yaml'))],

    # The __version__ comes from the exec at the top

    version=__version__,

    description="A blind search algorithm to look for transients in XMM and Chandra data",

    author='Giacomo Vianello',

    author_email='giacomo.vianello@gmail.com',

    url='https://github.com/giacomov/XtDac',

    download_url='https://github.com/giacomov/xdtac/archive/v%s' % __version__,

    keywords=['Blind search', 'Transients', 'Time-domain astronomy', 'X-rays'],

    classifiers=[],

    install_requires=[
        'numpy >= 1.6',
        'astropy >= 1.0',
        'scipy>=0.13',
        'fitsio',
        'pyyaml',
        'scikit-learn',
        'pyregion'],

    ext_modules=[],
    scripts = xmm_scripts + chandra_scripts + main_executables,

    package_data={
              'XtDac': ['data/chandra_csc_1.1.pickle.gz', 'data/sample_configuration.yml','data/obs_in_csc1.1.txt'],
           },
    include_package_data=True,
)
