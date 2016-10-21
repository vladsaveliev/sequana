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

Snakemake variant calling pipeline based on pipelines of Varun Khanna 
(https://github.com/khannavarun) and Adrien Villain (https://github.com/avillain).
Reads (paired or single) are mapped using bwa mem. Aligned reads are processed with
picard tools markduplicates and GATK indel realigner. Freebayes is used to detect
SNPs and short INDELs. An annotation file can be set to annotate detected variants.
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
