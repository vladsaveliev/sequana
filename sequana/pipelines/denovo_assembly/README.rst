:Overview: Denovo Assembly from FastQ files
:Input: FastQ file(s) from Illumina Sequencing instrument
:Output: FastA, VCF and HTML files


Usage
~~~~~~~~~

Example::

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

Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documenteted configuration file :download:`../sequana/pipelines/denovo_assembly/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file. Here are the rules and their developer and user documentation.

samtools_depth
^^^^^^^^^^^^^^
.. snakemakerule:: samtools_depth

bwa
^^^^
.. snakemakerule:: bwa_mem_dynamic

digital_normalisation
^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: digital_normalisation

format_contigs
^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: format_contigs

freebayes
^^^^^^^^^^
.. snakemakerule:: freebayes

quast
^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: quast

snpeff
^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: snpeff

spades
^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: spades
