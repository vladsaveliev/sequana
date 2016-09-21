:Overview: compressor can be used to compress/uncompress files recursively. It was designed to uncompress gzipped files and to compress them back into a bzip2 format. Supported formats are gz and bz2.
:Input: Any number of fastQ files (compressed or not)
:Output:
:requirements: gzip and bzip and their parallel version pigz and pbzip2

Usage
~~~~~~~

A standalone named sequana_compressor is provided. Please see::

    sequana_compressor --help 

Example::

    sequana_compressor --source fastq.gz --target fastq.bz2


Requirements
~~~~~~~~~~~~~~~~~~
Parallel version of gzip and bzip2::

.. include:: ../sequana/pipelines/compressor/requirements.txt


Details
~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/compressor/dag.png
