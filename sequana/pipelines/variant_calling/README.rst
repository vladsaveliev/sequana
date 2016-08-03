:Overview: Variant calling
:Input: fastq file from Illumina Sequencing instrument
:Output: vcf and html files
:Config file requirements:
    - samples:file1
    - samples:file2
    - project
    - reference:reference.fasta

Usage
~~~~~~~~~

::

    sequana --pipeline variant_calling --file1 R1.fastq.gz --file2 R2.fastq.gz --project variant --reference reference.fasta
    cd variant
    snakemake -s variant_calling.rules -p --stats stats.txt -j 4

Requirements
~~~~~~~~~~~~~~~~

- Bedtools
- BWA
- Freebayes
- GenomeAnalysisTK
- Picard-tools
- Samtools
- snpEff

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/variant_calling/variant_calling_dag.png
