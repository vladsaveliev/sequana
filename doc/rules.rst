
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


fastqc
==========

.. include:: ../sequana/rules/fastqc/README.rst


codecs
===========
gz_to_fastq
--------------
.. snakemakerule:: gz_to_fastq
bz2_to_fastq
-------------
.. snakemakerule:: bz2_to_fastq

.. snakemakerule:: dsrc_to_fastq
.. snakemakerule:: fastq_to_gz
.. snakemakerule:: fastq_to_bz2
.. snakemakerule:: fastq_to_dsrc


