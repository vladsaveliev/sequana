
.. _pipelines:

Pipelines
##############

This sections provides a summary of the official Sequana pipelines. Note that quality_taxon
re-uses quality, which itself uses phe phix_removal pipeline.

Each pipeline can be automatically downloaded (with its config file) using::

    sequana --pipeline <name> --file1 <FILE1> --file2 <FILE2>

.. contents::
    :depth: 2

.. _pipeline_phix_removal:



.. _denovo_assembly_pipeline:

denovo_assembly
====================

.. include:: ../sequana/pipelines/denovo_assembly/README.rst


Phix Removal
==============


.. include:: ../sequana/pipelines/phix_removal/README.rst

.. _pipeline_quality:

Quality
==========

.. include:: ../sequana/pipelines/quality/README.rst

.. _pipeline_quality_taxon:

Quality Taxon
=======================

.. include:: ../sequana/pipelines/quality_taxon/README.rst


.. _pipeline_variant_calling:

Variant Calling
=================

.. include:: ../sequana/pipelines/variant_calling/README.rst


.. _pipeline_compressor:

Compressor
============

.. include:: ../sequana/pipelines/compressor/README.rst
