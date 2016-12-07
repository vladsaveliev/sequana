Sequana documentation
##########################################

|version|, |today|


.. raw:: html

    <div style="width:80%"><p>    <a href="https://pypi.python.org/pypi/sequana"> <img src="https://badge.fury.io/py/sequana.svg"></a>
    <a href="https://travis-ci.org/sequana/sequana"> <img src="https://travis-ci.org/sequana/sequana.svg?branch=master"></a>
    <a href="https://coveralls.io/github/sequana/sequana?branch=master"> <img src="https://coveralls.io/repos/github/sequana/sequana/badge.svg?branch=master"></a>
    <a href="http://sequana.readthedocs.org/en/master/?badge=master"> <img src="http://readthedocs.org/projects/sequana/badge/?version=master"></a>
    </p><p>
    <a href="https://microbadger.com/images/sequana/sequana" title="Get your own image badge on microbadger.com">
    <img src="https://images.microbadger.com/badges/image/sequana/sequana.svg"></a>
    <a href="https://microbadger.com/images/sequana/sequana" title="Get your own version badge on microbadger.com"><img
    src="https://images.microbadger.com/badges/version/sequana/sequana.svg"></a>
    </p>
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
#. A set of :ref:`pipelines <Pipelines>` dedicated to NGS in the form of Snakefiles 
   (Makefile-like with Python syntax based on snakemake framework).
#. Original tools to help in the creation of such pipelines including HTML reports.
#. :ref:`Standalone applications<applications>`:
    #. :ref:`sequana_coverage<standalone_sequana_coverage>` ease the 
       extraction of genomic regions of interest and genome coverage information

Currently, the available pipelines cover quality control (e.g. adapters removal,
phix removal, trimming of bad quality bases), variant calling, characterisation
of the genome coverage, taxonomic classification, de-novo assembly. See the :ref:`pipelines`
section for more information.

**Sequana** can be used by developers to create new pipelines and by users in the
form of applications ready for production.


To join the project, please let us know on `github
<https://github.com/sequana/sequana/issues/306>`_.



.. Here we are building the carrousel

.. |bam| image::
      _images/sphx_glr_plot_bam_001.png
   :target: auto_examples/plot_bam.html

.. |coverage| image::
      _images/sphx_glr_plot_coverage_001.png
   :target: auto_examples/plot_coverage.html

.. |fastqc| image::
      _images/sphx_glr_plot_fastqc_hist_001.png
   :target: auto_examples/plot_fastqc_hist.html

.. |kraken| image::
      _images/sphx_glr_plot_kraken_001.png
   :target: auto_examples/plot_kraken.html


.. raw:: html

    <div class="body">
   <div id="index-grid" class="section group">
    <div class="col span_1_of_3">
        <h3><a href="installation.html">Installation</a></h3>
        <p>Using conda or docker</p>
        <h3><a href="auto_examples/index.html">Examples</a></h3>
        <p>Visit our example gallery to use the Python library</p>
        <h3><a href="pipelines.html">NGS pipelines</a></h3>
        <p>Learn about available Snakemake pipelines</p>
        <h3><a href="applications.html">Standalone applications</a></h3>
        <p>Standalone applications</p>
    </div>
    <div class="col span_2_of_3">
    <div class="jcarousel-wrapper">
    <div class="jcarousel">

* |coverage|
* |fastqc|
* |kraken|
* |bam|

.. raw:: html

            </div>
        <a href="#" class="jcarousel-control-prev">&lsaquo;</a>
        <a href="#" class="jcarousel-control-next">&rsaquo;</a>
        <p class="jcarousel-pagination">
        </p>
        </div>
        </div>
        </div>
   </div>
   <div style="clear: left"></div>









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

    sequana --pipeline quality_control --output-directory TEST --adapters PCRFree

This command downloads the required configuration file(s) in particular
the config file and the pipeline itself. This example should work out of
the box but you may want to look at the
configuration file **config.yaml**. For instance, you may want to change the
reference to the *phix* (by default we use *phix174.fa*, which is provided in Sequana) or
change the adapter_removal section to your needs (cutadapt parameters, in
particular the forward and reverse complement list of adapters; None by default).

By default, the output directory is called **analysis** and ca be overwritten with the ``--output-directory`` parameter. Then, run the pipeline and wait for completion.::

    cd TEST
    snakemake -s quality_control.rules --stats stats.txt -p -j 4 --forceall

The -p option shows the commands, -j 4 means use 4 threads when possible.
Alternatively, there is also a **runme.sh** script.

You should now have a directory with a HTML report correspinding to the sample::

    open Hm2_GTGAAA_L005/report_qc_Hm2_GTGAAA_L005/summary.html

More information can be found in the next sections.

User guide and reference
###########################


.. toctree::
    :numbered:
    :maxdepth: 2

    installation.rst
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
