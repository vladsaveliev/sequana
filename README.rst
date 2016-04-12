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

If you already already install dependencies, this should install the latest release::

    pip install sequana --upgrade

Some dependencies required include matplotlib, pandas, cutadapt, pysam. If you
are new or starting with Python, we strongly recommand you to use anaconda. We use the **bioconda** channel, which can be
added to your environment as follows::

    conda config --add channels r
    conda config --add channels bioconda
    
Then, install those dependencies::

    conda install numpy matplotlib pandas cutadapt pysam bwa bcftools pyvcf samtools snakemake biokit bioservices bedtools picard freebayes

although the code is Python2.7 and Python3.5 compatible, a dependency
(Snakemake) only supports Python3.5 for the moment so we will support only Python3.5 version. If you wish to use functionalities of Sequana that do not make use of Snakemake, you may still used it with Python2.7.

Developers can also install other tols::

    conda install nose coverage

**Please see the** `documentation <http://sequana.readthedocs.org>`_ **for usage and examples**





