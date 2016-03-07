# -*- coding: utf-8 -*-
__revision__ = "$Id: $" # for the SVN Id
import sys
import os
from setuptools import setup, find_packages
import glob

_MAJOR               = 0
_MINOR               = 0
_MICRO               = 2
version              = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
release              = '%d.%d' % (_MAJOR, _MINOR)

metainfo = {
    'authors': {"main": ("yourname", "email@whatever.org")},
    'version': version,
    'license' : 'GPL',
    'download_url' : ['http://pypi.python.org/pypi/sequana'],
    'url' : ["http://github.com/sequana/"],
    'description': "Put a short description here" ,
    'platforms' : ['Linux', 'Unix', 'MacOsX', 'Windows'],
    'keywords' : [''],
    'classifiers' : [
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics']
    }

# files in share/data
datadir = os.path.join('pipelines')
datafiles = [(datadir, [f for f in glob.glob(os.path.join(datadir, '*'))])]


setup(
    name             = "sequana",
    version          = version,
    maintainer       = metainfo['authors']['main'][0],
    maintainer_email = metainfo['authors']['main'][1],
    author           = metainfo['authors']['main'][0],
    author_email     = metainfo['authors']['main'][1],
    long_description = open("README.rst").read(),
    keywords         = metainfo['keywords'],
    description      = metainfo['description'],
    license          = metainfo['license'],
    platforms        = metainfo['platforms'],
    url              = metainfo['url'],
    download_url     = metainfo['download_url'],
    classifiers      = metainfo['classifiers'],

    # package installation
    packages = find_packages(exclude=["test"]),

    include_package_data=True,

    install_requires = ["easydev>=0.9.12", "reports", "matplotlib", "pandas",
        "cutadapt==1.9.1", "pysam"],

    # uncomment if you have share/data files
    # package_data = {
    #     '': ['Snakefile*']
    # },

    zip_safe=False,
    #use_2to3 = True, # causes issue with nosetests
    entry_points = {
        'console_scripts':[
            'fastq_head=scripts.fastq_head:main',
            'fastq_count=scripts.fastq_count:main',
        ]
    },

)
