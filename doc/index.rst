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
:Issues: Please fill a report on `github/sequana/sequana/issues <https://github.com/sequana/sequana/issues>`_


What is Sequana ?
=====================

Sequana is a versatile tool that provides (i) pipelines dedicated to NGS in the form of Snakefiles (Makefile-like with Python syntax) (ii) original tools to help in the creation of such pipelines (.e.g., plotting, statistical analysis of results) (iii) HTML reports and (iv) standalone applications.

The pipelines available cover quality control (e.g. removing adapters, phix, trimming of 
bad quality bases) , variant calling, characterisation of the genome coverage,
and taxonomic classification. See the :ref:`pipelines`  section for more
information.

Sequana can be used by developers to create new pipelines and by users in the
form of applications ready for production.


Would you like to join us, please let us know on the github website. 

.. _installation:

Installation
#################


If you have already installed **Sequana** dependencies, this command
should install the latest release posted on Pypi website::

    pip install sequana --upgrade

If not, be aware that Sequana relies on a few dependencies that needs 
to be compiled (time consumming and requires proper C compilator). 
For instance, we use Matplotlib, Pandas, cutadapt but some pipelines 
also require more specific tools (e.g. BWA for read alignment). We therefore
strongly recommend to use Anaconda and in particular the **bioconda** channel, which can be
added to your environment as follows (once Anaconda is installed)::

    conda config --add channels r
    conda config --add channels bioconda

Here is a non exhaustive list of dependencies that should be enough to run the
current pipelines (version 0.1.4). We split the command on several lines to
emphasize the standard Anaconda packages and te bioconda ones but you
can use only one::

    conda install numpy matplotlib pandas cutadapt pysam pyvcf 
    conda install snakemake biokit bioservices
    conda install bwa bcftools samtools bedtools picard freebayes fastqc


.. note:: Sequana is not fully compatible with Python 2.7 since a dependency
    (Snakemake) is only available for Python 3.5. However, many core
    functionalities would work under Python 2.7 


.. _quick_start:

Quick start example: the quality pipeline
#############################################

The following example will show how to run the quality control on a pair of FastQ files.
The data comes from a sequencing (using HiSeq technology) of a
Measles virus. For testing purposes, you can download :download:`R1
<../sequana/resources/datag/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory. 

First, run the sequana standalone application to initialise the pipeline
**quality**::

    sequana --init quality --input-dir .  --project TEST

This download required file(s) in particular the config file and the pipeline
itself. This example should work out of the box but you may want to look at the
configuration file **config.yaml**. For instance, you may want to change the
reference to the *phix* (by default we use Coliphage_phix174, which is provided in Sequana).

Then, run the pipeline and wait for completion. Check the output for warnings or
errors::

    snakemake -s Snakefile --stats stats.txt -p -j 4 --forceall

The -p option shows the commands, -j 4 means use 4 threads when possible. Finally, open 
the report in ./TEST::

    open TEST/report_quality.html

This will copy the Snakefile locally as well as the corresponding configuration
file.


User guide and reference
###########################


.. toctree::
    :numbered: 
    :maxdepth: 2

    userguide
    tutorial
    pipelines
    auto_examples/index
    developers
    references
    faqs
    Changelog
    glossary


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
