
.. _rules:

Rules
##############

Here below we present a small subset of rules used within the different
pipelines. Rules can be incorporated within a pipeline as follows (given its
name, here e.g. *data*)::

    from sequana import snaketools as sm
    include: sm.module['data']


.. contents::
    :depth: 3

data
==============

.. include:: ../sequana/rules/data/README.rst


kraken contaminant
=====================
.. snakemakerule:: kraken

fastqc
==========

.. include:: ../sequana/rules/fastqc/README.rst

codecs
===========
The following rules are used by the :ref:`pipeline_compressor` standalone
application.

gz_to_fastq
--------------
.. snakemakerule:: gz_to_fastq

bz2_to_fastq
-------------
.. snakemakerule:: bz2_to_fastq

dsrc_to_fastq
--------------

.. snakemakerule:: dsrc_to_fastq

fastq_to_gz
--------------
.. snakemakerule:: fastq_to_gz

fastq_to_bz2
-------------
.. snakemakerule:: fastq_to_bz2

fastq_to_dsrc
---------------
.. snakemakerule:: fastq_to_dsrc

gz_to_bz2
---------------
.. snakemakerule:: gz_to_bz2

bz2_to_gz
---------------
.. snakemakerule:: bz2_to_gz




bwa_mem
===========
.. snakemakerule:: bwa_mem_dynamic




