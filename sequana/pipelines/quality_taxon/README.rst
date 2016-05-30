:Overview: On top of quality pipeline, we add a taxonomy search
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: 
    - <PROJECT>/cutadapt/<PROJECT>R1.cutadapt.fastq.gz
    - <PROJECT>/cutadapt/<PROJECT>R2.cutadapt.fastq.gz
    - >PROJECT>/kraken/kraken_to_krona.html
:Config file requirements:
    - samples:file1
    - samples:file2
    - project:
    - bwa_mem:reference
    - kraken:database


Usage
~~~~~~~

::

    sequana init quality_taxon --file1 R1.fastq.gz --file2 R2.fastq.gz
    python sequana_quality.py


Requirements
~~~~~~~~~~~~~~~~~~

- bwa
- samtools
- kraken


.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/quality_taxon/dag.png
