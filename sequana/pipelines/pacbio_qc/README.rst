:Overview: Quality control for Pacbio bam data
:Input: A set of bam files
:Output: report_{sample}.html where sample is the name of an input file

Usage
~~~~~~~

First copy the BAM files into a directory. Start sequanix in that directory.
Input BAM files must have the extension **.bam**.

Then, start sequanix as follows::

    sequana --pipeline pacbio_qc --input-directory . --working-directory analysis

In the configuration panel, set the Kraken database. 

Save the project and press Run. Once done, open the HTML report for the bam of
interest



Requirements
~~~~~~~~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/pacbio_qc/dag.png


Details
~~~~~~~~~


This pipeline takes as inputs a set of BAM files from Pacbio sequencers. It
computes a set of basic statistics related to read length. It also shows some 
histogram related to the GC content, SNR of the diodes and the so-called ZMW
values.
