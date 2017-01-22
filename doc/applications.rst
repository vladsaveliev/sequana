
.. _applications:

Applications (standalone)
============================

.. contents::


.. _sequanix:

Sequanix
----------

:Overview: a Graphical User Interface (GUI) for Sequana pipelines and any
    Snakemake-based workflows.

Just type **sequanix** in a shell. 

.. note:: tested under Linux only. However, Mac and Windows users should be 
   able to use it since it is based on Python and PyQt. Again, we strongly
   advice to use Anaconda to install all required dependencies

Here is a snapshot.

.. image:: _static/sequanix.png


sequana
---------

:Overview: Creates project(s) to run a **Sequana** pipeline(s)

The **sequana** executable can be used to create pipelines (and associated
config file). For example::

    sequana --pipeline quality --file1 R1.fastq.gz --file2 R2.fastq.gz --project TEST

will create a directory called TEST with a few files such as *quality.rules*,
*config.yaml*, a *runme.sh* and a *README* file.

Valid pipelines can be found using::

    sequana --show-pipelines

There are many more options and documentation. Please use the ``--help``
option for more information.

.. _standalone_sequana_coverage:

sequana_coverage
--------------------

:Description: Show coverage and interval of confidence to identify under and
    over represented genomic regions.
:Help: please use sequana_coverage --help
:Docker: ::
    
        git pull sequana/sequana_coverage 

    See `github sequana_coverage docker page <https://github.com/sequana/sequana/tree/master/docker/sequana_coverage>`_ for details
:Sequana: See :class:`~sequana.bedtools.GenomeCov` to use the coverage in your own script.
:Gallery: See examples in the `gallery <http://sequana.readthedocs.io/en/master/auto_examples/index.html>`_

Starting from a BED file and its reference, one can use this command in a
shell::

    sequana_coverage  --input JB409847.sorted.bed -o
                      --reference JB409847.fa --show-html

It creates an HTML report with various images showing the coverage and GC
versus coverage plots. It also provides a set of CSV files with low or high
coverage regions (as compared to the average coverage).

.. seealso:: the underlying algorithm is described in details in the documentation
    (:mod:`sequana.bedtools.GenomeCov`).


sequana_summary
------------------

:Description: Prints basic statistics about a set of NGS input files. Currently
    handles Fastq (gzipped or not) or BED files (coverage).


sequana_mapping
------------------
:Description: a simple application to map reads onto a genome given one or two
    FastQ files (gzipped) and a refenrece.


sequana_taxonomy
--------------------

:Description: Creates a HTML document with Krona and pie chart of taxonomic
    content of a st of FastQ files. Uses Kraken and a dedicated Sequana
    database.

fastq related: fastq_count
-----------------------------

:Description: count number of reads and lines

fastq related: fastq_head 
-----------------------------

:Description: Extract head of a fastq files (zipped or not)


sequana_compressor
---------------------

:Description: converts fastq into a *gz* or *bz2* or *dscr* format. Conversely,
    can decompress or even convert from one compressed format to another one.
    Recursivit is allowed and thanks to snakemake pipeline used behind the scene, it
    works on a cluster and in parallel mode. See::

        sequana_compressor --help

