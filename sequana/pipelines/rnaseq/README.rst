:Overview: RNASeq
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: BAM, count and HTML files



Usage
~~~~~~~~~

Example::

    sequana --pipeline rnaseq --input-dir . --output-directory analysis --no-adapters
    cd analysis
    srun snakemake -s rnaseq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json --cluster "sbatch --mem={cluster.ram} --cpus-per-task={threads}"

Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/rnaseq/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/rnaseq/dag.png


Details
~~~~~~~~~

Snakemake RNA-seq pipeline based on workflow use at Biomics Pole in Institut Pasteur. The pipeline runs some QC, such as FastQC, fastq_screen (you need your own base).
Reads could be trimmed by several tools (cutadapt, atropos, clean_ngs) and mapped against reference genome (with bowtie or STAR) and ribosomal RNA (with bowtie1).
After, reads are counted with feature-counts (HTSeq-count soon available) against a GFF file. All results are summarized by multiQC.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/rnaseq/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. Here are the rules and their developer and user documentation.



FastQC
^^^^^^^^^^^
.. snakemakerule:: fastqc

Fastq_screen
^^^^^^^^^^^^^^^
.. snakemakerule:: fastq_screen

Cutadapt
^^^^^^^^^
.. snakemakerule:: cutadapt

Mapping on rRNA
^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: bowtie1_mapping_dynamic

Mapping on reference genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bowtie1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: bowtie1_mapping_dynamic

STAR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: star_mapping

Counting
^^^^^^^^^^^^
.. snakemakerule:: feature_counts

Reporting
^^^^^^^^^^^^
.. snakemakerule:: multiqc
