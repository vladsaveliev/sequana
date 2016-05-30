SEQUANA
############

.. image:: https://badge.fury.io/py/sequana.svg
    :target: https://pypi.python.org/pypi/sequana

.. image:: https://travis-ci.org/sequana/sequana.svg?branch=master
    :target: https://travis-ci.org/sequana/sequana

.. image:: https://coveralls.io/repos/github/sequana/sequana/badge.svg?branch=master
    :target: https://coveralls.io/github/sequana/sequana?branch=master 

.. image:: http://readthedocs.org/projects/sequana/badge/?version=master
    :target: http://sequana.readthedocs.org/en/latest/?badge=master
    :alt: Documentation Status

:Python version: Python 2.7 and 3.5
:Online documentation: `On readthedocs <http://sequana.readthedocs.org/>`_
:Issues and bug reports: `On github <https://github.com/sequana/sequana/issues>`_






**Sequana** includes a set of pipelines related to NGS (new generation sequencing) including quality control, variant calling, coverage, taxonomy. 

**Please see the** `documentation <http://sequana.readthedocs.org>`_ **for usage and examples**

Installation
=================

The installation process is explained in the documentation but here is a quick
explanation. If you have already installed **Sequana** dependencies, this command
should install the latest release posted on Pypi website::

    pip install sequana --upgrade

There are a few dependencies that needs to be compiled (time consumming and
requires proper C compilator). For instance, we use Matplotlib, Pandas, cutadapt but some pipelines also require more specific tools (e.g. BWA for read alignment). We therefore strongly recommend to use Anaconda and in particular the **bioconda** channel, which can be
added to your environment as follows (once Anaconda is installed)::

    conda config --add channels r
    conda config --add channels bioconda

Here is a non exhaustive list of dependencies that should be enough to run the
current pipelines (version 0.1.4). We split the command on several lines to
emphasize the standard Anaconda packages and te bioconda ones but you
can use only one::

    conda install numpy matplotlib pandas cutadapt pysam pyvcf snakemake biokit bioservices
    conda install bwa bcftools samtools bedtools picard freebayes fastqc


.. note:: Sequana is compatible with Python 2.7 and 3.5 but pipelines are
    built with Snakemake, which is only 3.5 compatible. So officially, Sequana is a
    Python3.5 compatible package only.







