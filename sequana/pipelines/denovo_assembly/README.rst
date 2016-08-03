:Overview: Denovo Assembly
:Input: fastq file from Illumina Sequencing instrument
:Output: fasta, bam, vcf and html files
:Config file requirements:
    - samples:file1
    - samples:file2
    - project

Usage
~~~~~~~~~

::

    sequana --pipeline denovo_assembly --file1 R1.fastq.gz --file2 R2.fastq.gz --project denovo
    cd denovo
    snakemake -s denovo_assembly.rules -p --stats stats.txt -j 4
    

Requirements
~~~~~~~~~~~~~~~~

- bedtools
- bwa
- freebayes
- genomeAnalysisTK
- khmer
- picard-tools
- quast
- samtools
- snpEff
- spades

Details
~~~~~~~~~

The reads normalization with khmer (digital normalization) 
is an optional task setting up in the config file.
You might use normalization if the depth of coverage is too important.
Then, the denovo assembly is done with spades.
Finally, the coverage and misassembly are evaluated with the variant calling
pipeline.


.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/denovo_assembly/denovo_dag.png
