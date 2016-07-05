# -*- coding: utf-8 -*-
# License: 3-clause BSD
__revision__ = "$Id: $" # for the SVN Id
import sys
import os
from setuptools import setup, find_packages
import glob

_MAJOR               = 0
_MINOR               = 1
_MICRO               = 6
version              = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
release              = '%d.%d' % (_MAJOR, _MINOR)

metainfo = {
    'authors': {"main": ("yourname", "email@whatever.org")},
    'version': version,
    'license' : 'new BSD',
    'download_url' : ['http://pypi.python.org/pypi/sequana'],
    'url' : ["http://github.com/sequana/"],
    'description': "Put a short description here" ,
    'platforms' : ['Linux', 'Unix', 'MacOsX', 'Windows'],
    'keywords' : [''],
    'classifiers' : [
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics']
    }


packages = find_packages()
packages = [this for this in packages if this.startswith('test.') is False]
packages = [this for this in packages if this not in ['test']]

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
    packages = packages,

    # pillow, sphinx-gallery and numpydoc are  for the doc only
    # mock is for the test only
    install_requires = ["easydev>=0.9.22", "reports>=0.1.8", "matplotlib",
        "pandas", "cutadapt>=1.9.1", "bioservices>=1.4.11", "biokit>=0.3.0", 
        "pysam", "pyVCF", "sphinx-gallery", "mock", "numpydoc", "pillow"],

    # here below '': pattern means include that pattern in all packages
    # so '' :['README.rst'] will include all README.rst recursively
    # required to use python setup.py install

    # This is recursive include of data files
    exclude_package_data = {"": ["__pycache__"]},
    package_data = {
        '': ['Snakefile*', '*html', 'README.rst', 'config.yaml*', '*.css', "*.js", 
                "snpEff.config*", "*.fa", "*.rules"],
        'sequana.rules' : ['*/*.rules', "*/*/*.rules"],
        'sequana.pipelines' : ['*/*.rules', "*/*/*.rules"],
        'sequana.resources.data' : ['*'],
        'sequana.resources.images' : ['*'],
        'sequana.resources.testing' : ['*'],
        'sequana.resources.js/galleria/themes' : ['*'],
        },

    # thise files do not need to be added in MANIFEST.in since there are python
    # packages that will be copied from sequana/ into sequana/
    # Note, however, that e.g. ./pipelines must be added 

    zip_safe=False,
    entry_points = {
        'console_scripts':[
           'fastq_head=sequana.scripts.fastq_head:main',
           'fastq_count=sequana.scripts.fastq_count:main',
           'sequana=sequana.scripts.main:main',
           'sequana_taxonomy=sequana.scripts.taxonomy:main',
           'sequana_coverage=sequana.scripts.coverage:main'
        ]
    },

)
