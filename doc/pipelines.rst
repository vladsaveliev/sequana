
.. _pipelines:

Pipelines
##############

This sections provides a summary of the official Sequana pipelines, which are all based on snakemake.

Each pipeline can be automatically downloaded (with its config file) using::

    sequana --pipeline <name> --file1 <FILE1> --file2 <FILE2>

.. contents::
    :depth: 2


.. _denovo_assembly_pipeline:

denovo_assembly
====================

.. include:: ../sequana/pipelines/denovo_assembly/README.rst


.. _pipeline_quality_control:

Quality control
=================

.. include:: ../sequana/pipelines/quality_control/README.rst

.. _pipeline_rnaseq:

RNA-Seq 
=================

.. include:: ../sequana/pipelines/rnaseq/README.rst


.. _pipeline_variant_calling:

Variant Calling
=================

.. include:: ../sequana/pipelines/variant_calling/README.rst


.. _pipeline_compressor:

Compressor
============

.. include:: ../sequana/pipelines/compressor/README.rst
