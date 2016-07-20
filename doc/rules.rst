
.. _rules:

Rules
##############

Here below we present a small subset of rules used within the different
pipelines. Rules can be incorporated within a pipeline as follows (given its
name, here e.g. *data*)::

    from sequana import snaketools as sm
    include: sm.module['data']


.. contents::
    :depth: 2

data
==============

.. include:: ../sequana/rules/data/README.rst


fastqc
==========

.. include:: ../sequana/rules/fastqc/README.rst




