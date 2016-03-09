SEQUANA
############

.. image:: https://badge.fury.io/py/sequana.svg
    :target: https://pypi.python.org/pypi/sequana

.. image:: https://travis-ci.org/sequana/sequana.svg?branch=master
    :target: https://travis-ci.org/sequana/sequana

.. image:: https://coveralls.io/repos/github/sequana/sequana/badge.svg?branch=master
    :target: https://coveralls.io/github/sequana/sequana?branch=master 

.. image:: http://readthedocs.org/projects/sequana/badge/?version=latest
    :target: http://sequana.readthedocs.org/en/latest/?badge=latest
    :alt: Documentation Status

:Python version: Python 2.7, 3.4 and 3.5
:Online documentation: `On readthedocs <http://sequana.readthedocs.org/>`_
:Issues and bug reports: `On github <https://github.com/sequana/sequana/issues>`_






**Sequana** includes a set of pipelines related to NGS (new generation sequencing). 

It will provide a set of modular pipelines and reports associated to them.


Installation
=================


::

    pip install sequana


Some dependencies required include matplotlib, pandas, cutadapt, pysam. If you
are new or starting with Python, we strongly recommand to use anaconda and to
install those dependencies::

    conda install matplotlib pandas cutadapt pysam

although the code is Python2.7 and Python3.5 compatible, a dependency
(Snakemake) only supports Python3.5 for the moment so, we will assume you have a
Python3.5 version installed.


Example
==========

::

    # contaminant and adapter removal search
    from sequana import biomics


    pipeline = biomics.Biomics("*fastq.gz") # assumm two paired-end fastq files
    pipeline.create()
    pipeline.start()
    pipeline.report()

Behind the scene, a Snakefile and a config file are created in the local
directory. The snakefile can be executed indepedently::

    snakemake Snakefile -p




