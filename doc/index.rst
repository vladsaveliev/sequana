Sequana documentation
##########################################

|version|, |today|


.. image:: https://badge.fury.io/py/sequana.svg
    :target: https://pypi.python.org/pypi/sequana

.. image:: https://travis-ci.org/sequana/sequana.svg?branch=master
    :target: https://travis-ci.org/sequana/sequana

.. image::  https://coveralls.io/repos/github/sequana/sequana/badge.svg?branch=master
   :target: https://coveralls.io/github/sequana/sequana?branch=master 

.. image:: http://readthedocs.org/projects/sequana/badge/?version=master
    :target: http://sequana.readthedocs.org/en/master/?badge=master
    :alt: Documentation Status


:Python version: Python 3.5 although some modules are Python2.7 compatible
:Issues: Please fill report on `github/sequana/sequana/issues <https://github.com/sequana/sequana/issues>`_


What is Sequana ?
=====================

Sequana is a versatile tool that provides (i) pipelines dedicated to NGS in the form of Snakefiles (Makefile-like with Python syntax) (ii) original tools to help in the creation of such pipelines (.e.g., plotting, statistical analysis of results) (iii) HTML reports and (iv) standalone application.

Sequana can be used by developers to create new pipelines and by users in the
form of applications ready for production.


Would you like to join us, please let us now on the github website. 


Installation
#################

If you have a Python environment, you may use:: 

    pip install sequana

Otherwise, we would recommend to use Conda and install these dependencies
first::

    conda install numpy pandas matplotlib
    pip install sequana

Some pipelines have specific dependencies that can be installed with the
**bioconda** channel.

Examples
##########

Given a set of FASTQ files in a local directory, you can run a pipeline as
follows (given its name)::

    sequana --init quality --file1 R1.fastq.gz --file2 R2.fastq.gz --project TEST
    # EDIT the file config.yaml to provide further information such as the
    # kraken database or possibly reference for a phix (by default uses
    # Coliphase_phix174

Run the pipelines and wait for completion. Check the output ::

    snakemake -s Snakefile --stats stats.txt -p -j 4 --forceall

The -p option shows the commands, -j 4 means use 4 threads when possible.

Finally, open the report in ./TEST::

    open TEST/report_quality.html


This will copy the Snakefile locally as well as the corresponding configuration
file.


User guide and reference
###########################


.. toctree::
    :numbered: 
    :maxdepth: 2

    userguide
    pipelines
    auto_examples/index
    references
    faqs
    Changelog


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
