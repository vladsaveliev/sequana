:Overview: Quality control for Pacbio bam data (raw data)
:Input: A set of bam files
:Output: report_{sample}.html where sample is the name of an input file

Usage
~~~~~~~

First copy the BAM files into a directory. Start **Sequanix** in that directory.
Input BAM files must have the extension **.bam**.

Then, start sequanix as follows::

    sequanix -w analysis -i . -p pacbio_qc

In the sequana/input file panel, set the extension to `*.bam`

In the configuration tab, in the kraken section add as many databases 
as you wish.

Save the project and press Run. Once done, open the HTML report for the BAM of
interest.


Requirements
~~~~~~~~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/pacbio_qc/dag.png


Details
~~~~~~~~~

This pipeline takes as inputs a set of BAM files from Pacbio sequencers. It
computes a set of basic statistics related to read length. It also shows some 
histogram related to the GC content, SNR of the diodes and the so-called ZMW
values. Finally, a quick taxonomy can be performed.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/pacbio_qc/config.yaml` to be used with the pipeline.

.. snakemakerule:: bam_to_fasta
