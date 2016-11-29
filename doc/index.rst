Sequana documentation
##########################################

|version|, |today|


.. raw:: html

    <div style="width:80%"><p>    <a href="https://pypi.python.org/pypi/sequana"> <img src="https://badge.fury.io/py/sequana.svg"></a> 
    <a href="https://travis-ci.org/sequana/sequana"> <img src="https://travis-ci.org/sequana/sequana.svg?branch=master"></a>
    <a href="https://coveralls.io/github/sequana/sequana?branch=master"> <img src="https://coveralls.io/repos/github/sequana/sequana/badge.svg?branch=master"></a>
    <a href="http://sequana.readthedocs.org/en/master/?badge=master"> <img src="http://readthedocs.org/projects/sequana/badge/?version=master"></a> 
    </p><p>
    <a href="https://microbadger.com/images/cokelaer/sequana"> <img src="https://images.microbadger.com/badges/image/cokelaer/sequana.svg"> </a> 
    <a href="https://microbadger.com/images/cokelaer/sequana"> <img src="https://images.microbadger.com/badges/version/cokelaer/sequana.svg"> </a></p>
    </div>


:Python version: Python 3.5 although some modules are Python2.7 compatible
:Source: See  `http://github.com/sequana/sequana <https://github.com/sequana/sequana/>`_.
:Issues: Please fill a report on `github <https://github.com/sequana/sequana/issues>`_
:How to cite: Desvillechabrol D, Bouchier C, Cokelaer T and Kennedy S. Sequana: a set of
    flexible genomic pipelines for processing and reporting NGS analysis [v1; no peer reviewed]. 
    `F1000Research 2016, 5:1767 <http://f1000research.com/posters/5-1767>`_ (poster) (doi:
    10.7490/f1000research.1112656.1)


What is Sequana ?
=====================

**Sequana** is a versatile tool that provides

#. A Python library dedicated to NGS analysis (e.g., tools to visualise standard NGS formats).
#. A set of :ref:`pipelines <Pipelines>` dedicated to NGS in the form of Snakefiles (Makefile-like with Python syntax based on snakemake framework).
#. Original tools to help in the creation of such pipelines including HTML reports.
#. :ref:`Standalone applications.<applications>`

Currently, the available pipelines cover quality control (e.g. adapters removal, 
phix removal, trimming of bad quality bases), variant calling, characterisation 
of the genome coverage, taxonomic classification, de-novo assembly. See the :ref:`pipelines`
section for more information.

**Sequana** can be used by developers to create new pipelines and by users in the
form of applications ready for production.


Would you like to join us, please let us know on the github website. 

.. _installation:

Installation
#################

If you have already installed **Sequana** dependencies, this command
should install the latest release posted on http://pypi.python.org/pypi/sequana website::

    pip install sequana --upgrade

If not, be aware that Sequana relies on many dependencies that needs
to be compiled (i.e., it is time consumming and requires proper C compilator).
For example, we use Matplotlib, Pandas, cutadapt that are Python libraries. 
However, many pipelines rely on third-party software such as BWA, Spades that are not
Python libraries. In practice, we do use `Anaconda <https://anaconda.org/>`_ and in
particular the **bioconda** channel, which can be
added to your environment as follows (once Anaconda is installed)::

    conda config --add channels r
    conda config --add channels bioconda

Some packages musts be installed::

    conda install numpy matplotlib pandas snakemake graphviz scipy

Then, depending on the pipelines or standalone applications you want to use,
you will need to install other packages. Here is a non exhaustive list
of dependencies that should be enough to run most of the current
pipelines (commands are split on several lines but you can also
install everything in one go)::

    conda install pysam pyvcf snpeff biokit bioservices spades khmer
    conda install bwa bcftools samtools bedtools picard freebayes fastqc
    conda install kraken krona pigz


.. note:: the denovo_assembly pipelines uses Quast tool. We ported the package to
   python 3.5 but this is not yet in bioconda. One can
   install it from the quast github (sept 2016).

.. note:: **Sequana** is not fully compatible with Python 2.7 since a dependency
    (Snakemake) is only available for Python 3.5. However, many core
    functionalities would work under Python 2.7

.. note:: We also provide a docker file with release >=0.14 of Sequana pre-installed.
    Please see https://github.com/sequana/sequana/tree/master/docker for
    details.

.. _quick_start:

Quick start example: the quality pipeline
#############################################

**Sequana** comes with standalone applications and pipelines in the form of
Snakefile (`snakemake <https://bitbucket.org/snakemake/snakemake/wiki/Home>`_)

The following example will show how to initialise and run the quality control pipeline
on a pair of FastQ files. 
The data comes from a sequencing (using HiSeq technology) of a
Measles virus. For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory. 

First, run the sequana standalone application to initialise the pipeline
**quality_control**::

    sequana --pipeline quality_control --input-dir . --output-directory TEST --adapters PCRFree

This command downloads the required configuration file(s) in particular 
the config file and the pipeline itself. This example should work out of 
the box but you may want to look at the
configuration file **config.yaml**. For instance, you may want to change the
reference to the *phix* (by default we use *phix174.fa*, which is provided in Sequana) or
change the adapter_removal section to your needs (cutadapt parameters, in
particular the forward and reverse complement list of adapters; None by default).

Note that the ``--output-directory`` parameter is optional. If not provided, the directory is named  **analysis**. 

Then, run the pipeline and wait for completion.::

    cd TEST
    snakemake -s quality_control.rules --stats stats.txt -p -j 4 --forceall

The -p option shows the commands, -j 4 means use 4 threads when possible.
Alternatively, there is also a **runme.sh** script.

Check the output and if there is no errors, you should be able to open the HTML report 
corresponding to the sample::

    open Hm2_GTGAAA_L005/report_qc_Hm2_GTGAAA_L005/summary.html

More information can be found in the next sections.

User guide and reference
###########################


.. toctree::
    :numbered: 
    :maxdepth: 2

    userguide
    tutorial
    pipelines
    auto_examples/index
    case_examples
    developers
    rules
    applications
    references
    faqs
    Changelog
    glossary


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
