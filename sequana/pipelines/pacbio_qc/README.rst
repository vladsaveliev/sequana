:Overview: Quality control for Pacbio bam data (raw data)
:Input: A set of bam files
:Output: report_{sample}.html where sample is the name of an input file

Usage
~~~~~~~

First copy the BAM files into a directory. Start **Sequanix** in that directory
as follows (input BAM files must have the extension **.bam**)::

    sequanix -w analysis -i . -p pacbio_qc

You are ready to go. If you want to filter out some BAM files, you may use the
pattern in tab 'input data'.

In the configuration tab, in the kraken section add as many databases
as you wish. You may simply unset the first database to skip the taxonomy, which
is experimental.

Save the project and press Run. Once done, open the HTML report for the BAM of
interest.


Requirements
~~~~~~~~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/pacbio_qc/dag.png


Details
~~~~~~~~~

This pipeline takes as inputs a set of BAM files from Pacbio sequencers. It
computes a set of basic statistics related to the read lengths. It also shows some
histograms related to the GC content, SNR of the diodes and the so-called ZMW
values. Finally, a quick taxonomy can be performed using Kraken. HTML reports
are created for each sample as well as a multiqc summary page.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/pacbio_qc/config.yaml` to be used with the pipeline.

bam_to_fasta
^^^^^^^^^^^^

.. snakemakerule:: bam_to_fasta

pacbio_quality
^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: pacbio_quality
