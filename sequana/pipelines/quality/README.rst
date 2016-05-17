:Overview: On top of phix_removal pipeline, remove adapters
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: 
    - <PROJECT>/cutadapt/<PROJECT>R1.cutadapt.fastq.gz
    - <PROJECT>/cutadapt/<PROJECT>R2.cutadapt.fastq.gz
:Config file requirements:
    - samples:file1
    - samples:file2
    - project:
    - bwa_mem:reference


Usage
~~~~~~~

::

    sequana init quality --file1 R1.fastq.gz --file2 R2.fastq.gz
    python sequana_quality.py


Requirements
~~~~~~~~~~~~~~~~~~

- bwa
- samtools



.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/quality/dag.png
