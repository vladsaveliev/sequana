:Overview: cleanup FastQ raw data from a Phix (known contaminant) using BWA
:Input: FastQ raw data (paired or single-end)
:Output: FastQ raw with clean reads (no phix) in ./<PROJECT>/bwa_bam_to_fastq/
:Config file requirements:    
    - samples:file1
    - samples:file2
    - project:
    - bwa_ref:reference


Usage
~~~~~~~

Example::

    sequana --pipeline phix_removal --glob "*gz" --project Phix
    cd Phix
    snakemake -s phix_removal -p --stats stats.txt

.. note:: Change the config.yaml file if needed 

Requirements
~~~~~~~~~~~~~~~~~~

- bwa
- samtools



.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/phix_removal/dag.png
