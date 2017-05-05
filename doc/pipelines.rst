
.. _pipelines:

Pipelines
##############

In **Sequana** parlance, a pipeline is an application based on Snakemake that consists of a Snakefile and a configuration file. Each pipeline and its configuration file can be automatically downloaded using::

    sequana --pipeline <name> 

By default the previous command creates a directory named *analysis* where the
pipeline and config file are stored. The pipeline must not be changed but the
configuration file can be edited to change the options. 

Although the configuration is documented and should be self content, additional
help for users and developers can be found for each pipeline in the following
links.


.. todo:: the following sections are in progress but should already give useful
   information about the pipelines that are available.

.. toctree::
    :maxdepth: 1

    pipeline_denovo_assembly.rst
    pipeline_quality_control.rst
    pipeline_rnaseq.rst
    pipeline_smallrnaseq.rst
    pipeline_variant_calling.rst
    pipeline_compressor.rst
    pipeline_pacbio_qc.rst


