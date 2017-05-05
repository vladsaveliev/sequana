:Overview: Quality control for Pacbio bam data
:Input: A set of bam files
:Output: report_{sample}.html where sample is the name of an input file

Usage
~~~~~~~

Please use sequanix.

::

    sequana --pipeline pacbio_qc --input-directory . --working-directory analysis --extension bam


Requirements
~~~~~~~~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/pacbio_qc/dag.png


Details
~~~~~~~~~


This pipeline takes as inputs a set of BAM files from Pacbio sequencers. It
computes a set of basic statistics related to read length. It also shows some 
histogram related to the GC content, SNR of the diodes and the so-called ZMW
values.
