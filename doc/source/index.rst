Sequana documentation
##########################################

|version|, |today|


.. image:: https://badge.fury.io/py/sequana.svg
    :target: https://pypi.python.org/pypi/sequana

.. image:: https://travis-ci.org/sequana/sequana.svg?branch=master
    :target: https://travis-ci.org/sequana/sequana

.. image::  https://coveralls.io/repos/github/sequana/sequana/badge.svg?branch=master
   :target: https://coveralls.io/github/sequana/sequana?branch=master 

.. image:: http://readthedocs.org/projects/sequana/badge/?version=latest
    :target: http://sequana.readthedocs.org/en/latest/?badge=latest
    :alt: Documentation Status

What is Sequana ?
=====================

A Python package that provides NGS pipelines in the form of Snakefiles (link to
snakemake). Those pipelines uses utilities and reports provided within Sequana
so as to ease the elaoration of new pipelines. Pipelines are made of
modules that we try to make generic enough to be re-used if possible.

Installation
#################


:: 

    pip install sequana

With conda::

    conda install numpy pandas matplotlib
    pip install sequana

Examples
##########

Given a set of FASTQ files in a local directory, you can run a pipeline as
follows (given its name)::

    sequana --pipeline contaminant

This will copy the Snakefile locally as well as the corresponding configuration
file.

Pipelines can be found in the directory ./pipelines



Issues
##########

Please fill bug report in https://github.com/sequana/sequana



.. toctree::
    :maxdepth: 2

    userguide
    references


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
