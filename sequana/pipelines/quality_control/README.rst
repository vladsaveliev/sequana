:Overview: Quality control, trimming (adapter removal) and taxonomic overview
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: 
    - <SAMPLE>/<SAMPLE>_report_qc/<SAMPLE>_R1_.cutadapt.fastq.gz
    - <SAMPLE>/<SAMPLE>_report_qc/<SAMPLE>_R1_.cutadapt.fastq.gz
    - <SAMPLE>/<SAMPLE>_report_qc/kraken/kraken.html
:Config file requirements:
    - samples: file1
    - samples: file2
    - project:
    - bwa_mem: reference
    - kraken: database


Usage
~~~~~~~

::

    sequana --pipeline quality_control --file1 R1.fastq.gz --file2 R2.fastq.gz
    python sequana_quality.py


Requirements
~~~~~~~~~~~~~~~~~~

- bwa
- samtools
- kraken
- krona

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/quality_control/dag.png
