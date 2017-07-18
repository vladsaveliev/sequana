:Overview: Variant calling
:Input: fastq file from Illumina Sequencing instrument
:Output: vcf and html files
:Config file requirements:
    - samples:file1
    - samples:file2
    - project
    - reference:reference.fasta

Details
~~~~~~~~

Snakemake variant calling pipeline based on
`tutorial <https://github.com/ekg/alignment-and-variant-calling-tutorial>`_
written by Erik Garrison. Reads (paired or single) are mapped using
`bwa <http://bio-bwa.sourceforge.net/>`_ and sorted with
`sambamba-sort <http://lomereiter.github.io/sambamba/docs/sambamba-sort.html>`_.
PCR duplicates are marked with
`sambamba-markdup <http://lomereiter.github.io/sambamba/docs/sambamba-sort.html>`_. 
`Freebayes <https://github.com/ekg/freebayes>`_ is used to detect SNPs and short
INDELs. Good practice proposed by the Broad institute is not necessary using
Freebayes as variant caller. Thus, the use of Genome Analysis Toolkit (GATK) is
avoided. For more information, please refer to a post by Brad Chapman

An annotation file can be set to annotate detected variants.
Variants are reported in a HTML report.

The pipeline provides a coverage analysis of the mapping coverage after bam processing.
Coverage for each base position is computed with bedtools genomecov. Sequana provides
a HTML report with dynamics plots of sequencing coverage and shows interesting regions 
which have unusual coverage depth.

Usage
~~~~~~~~~

::

    sequana --pipeline variant_calling --file1 R1.fastq.gz --file2 R2.fastq.gz --project variant --reference reference.fasta
    cd variant
    snakemake -s variant_calling.rules -p --stats stats.txt -j 4

Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/variant_calling/requirements.txt


.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/variant_calling/variant_calling_dag.png
