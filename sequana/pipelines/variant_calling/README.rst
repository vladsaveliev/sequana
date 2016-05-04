:Overview: Variant calling
:Input: fastq file from Illumina Sequencing instrument
:Output: bam, vcf and html files
:Config file requirements:
    - samples:file1
    - samples:file2
    - project
    - reference:sequence.fa

Usage
~~~~~~~~~

::

    mkdir analysis
    cd analysis
    # Put all necessary file in directory
    cp path/to/config.yaml .
    # Edit the config.yaml file
    snakemake -p -s path/to/Snakefile

Requirements
~~~~~~~~~~~~~~~~

- Bedtools
- BWA
- Freebayes
- GenomeAnalysisTK
- Picard-tools
- Samtools
