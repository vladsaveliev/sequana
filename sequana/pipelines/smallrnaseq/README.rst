:Overview: Small RNASeq (miRNA) analysis
:Input: FastQ raw data from Illumina Sequencer (single end only)
:Output: BAM, count and HTML files




Usage
~~~~~~~~~

Example::

    sequana --pipeline smallrnaseq --input-dir .  --output-directory analysis --adapters TruSeq
    cd analysis
    srun snakemake -s smallrnaseq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json --cluster "sbatch --mem={cluster.ram} --cpus-per-task={threads}" --restart-times 2


Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/smallrnaseq/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/smallrnaseq/dag.png


Details
~~~~~~~~~

This pipeline allows to map and count reads on mature and hairpin sequences (to download from miRBase) and perform some QC on data. All results are summarized by multiQC.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/smallrnaseq/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. Here are the rules and their developer and user documentation.



FastQC
^^^^^^^^^^^
.. snakemakerule:: fastqc_dynamic

Fastq_screen
^^^^^^^^^^^^^^^
.. snakemakerule:: fastq_screen

Cutadapt
^^^^^^^^^
.. snakemakerule:: cutadapt

Bowtie1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mapping on hairpin and mature sequences. The corresponding fasta files must be in the config file.

.. snakemakerule:: bowtie1_mapping_dynamic

Counting
^^^^^^^^^^^^
.. snakemakerule:: miRNA_count_dynamic

Reporting
^^^^^^^^^^^^
.. snakemakerule:: multiqc
