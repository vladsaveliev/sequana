:Overview: **compressor** can be used to compress/uncompress FastQ files recursively. 
   It was designed to uncompress gzipped files and to compress them back into a 
   bzip2 format. Supported formats are gz and bz2 and 
   dsrc (http://sun.aei.polsl.pl/dsrc/download.html).
:Input: Any number of fastQ files (compressed or not)
:Output: The input files (compressed or not)
:requirements: gzip and bzip and their parallel versions (pigz and pbzip2) as
    well as dsrc (DNA compression tool).

Usage
~~~~~~~

A standalone named **sequana_compressor** is provided. Please see::

    sequana_compressor --help

Example::

    sequana_compressor --source fastq.gz --target fastq.bz2


.. image:: _static/compressor_codecs.png
   :width: 60%

**compressor** allows one to go from one format in (fastq, fastq.gz, fastq.bz2,
fastq.dsrc) to any format in the same list. So this is a fully connected
network:

.. warning:: currently, this works only for fastq files and the **fastq** string
   must be provided before the extension

.. warning:: during the conversion, a .snakemake is created in each processed
   direcrory. If you interrupt the process, snakemake locks the directory. If
   you get an error message about locked directories, relaunch your previous
   command with --unlock to unlock the directories and start again

Requirements
~~~~~~~~~~~~~~~~~~
Parallel version of gzip and bzip, as well as dsrc:

.. include:: ../sequana/pipelines/compressor/requirements.txt

Config file
~~~~~~~~~~~~~~

In principle you should not use any config file if you use the standalone.
Note, however, the format of the underlying config file (for pbzip2 and pigz,
the number of threads is automatically set to the number of available threads).

::

    compressor:
        source: fastq.gz
        target: fastq.bz2
        threads: 4
        recursive: True
        verbose: True


DAG
~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/compressor/dag.png




Rules used by the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on the value of the target and source, only one rule is included in the
pipeline. For example if your source is **fastq.gz** and the target **fastq.bz2**, the
**gz_to_bz2** rule is included. Its documentation is here below:

.. snakemakerule:: gz_to_bz2

Others similar rules that convert from one compressed format to another
compressed formats are:

.. snakemakerule:: gz_to_dsrc
.. snakemakerule:: bz2_to_dsrc
.. snakemakerule:: bz2_to_gz
.. snakemakerule:: dsrc_to_gz
.. snakemakerule:: dsrc_to_bz2


