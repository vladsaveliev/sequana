:Overview: cleanup FastQ raw data from a Phix (known contaminant)
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: FastQ raw with clean reads (no phix) in ./<PROJECT>/bwa_bam_to_fastq/
:Config file requirements:
    - samples:file1
    - samples:file2
    - project
    - bwa_mem:reference


Usage
~~~~~~~

::

    mkdir analysis
    cd analysis
    # EDIT the config.yaml file 
    snakemake -p --stats stats.txt

Requirements
~~~~~~~~~~~~~~~~~~

- bwa
- samtools



.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/phix_removal/dag.png
