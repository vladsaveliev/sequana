:Overview: Variant calling
:Input: fastq file from Illumina Sequencer
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
