:Overview: Denovo Assembly from FastQ files
:Input: FastQ file(s) from Illumina Sequencing instrument
:Output: FastA, VCF and HTML files
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

.. include:: ../sequana/pipelines/denovo_assembly/requirements.txt

Details
~~~~~~~~~

The reads normalization is performed with khmer (digital normalization). It
is an optional task defined in the config file.
You might use normalization if the sequencing depth (depth of coverage) is too important.
Then, the denovo assembly is done with spades.
Finally, the coverage and misassembly are evaluated with the variant calling
pipeline.


.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/denovo_assembly/denovo_dag.png
