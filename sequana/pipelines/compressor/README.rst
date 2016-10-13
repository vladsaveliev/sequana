:Overview: compressor can be used to compress/uncompress files recursively. It was designed to uncompress gzipped files and to compress them back into a bzip2 format. Supported formats are gz and bz2 and dsrc (http://sun.aei.polsl.pl/dsrc/download.html).
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


Details
~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/compressor/dag.png
