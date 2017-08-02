:Overview: RNASeq: Differential expressed genes analysis
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

Snakemake RNA-seq pipeline based on workflow use at Biomics Pole in Institut Pasteur. The pipeline runs some QC, such as FastQC, fastq_screen (you need your own base).Reads could be trimmed by several tools (cutadapt, atropos, clean_ngs) and mapped against reference genome (with bowtie or STAR, bowtie2 is used by fastq_screen) and ribosomal RNA (with bowtie1). After, reads are counted with feature-counts (HTSeq-count soon available) against a GFF file. All results are summarized by multiQC.

.. warning:: The statistical analysis is not included in our pipeline because it is a step that is difficult to automate before to explore the data. You can continue with SARTools (https://github.com/PF2-pasteur-fr/SARTools) if you want to perform the statistical step.

Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/rnaseq/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. Here are the rules and their developer and user documentation.



FastQC
^^^^^^^^^^^

FastQC is used to check quality of sequenced reads.

.. snakemakerule:: fastqc

Fastq_screen
^^^^^^^^^^^^^^^

Fastq_screen is used to search any contamination in data. A interne database (with bowtie2 indexes) is mandatory.

.. snakemakerule:: fastq_screen

Cutadapt
^^^^^^^^^

Cutadapt is used to trim and filter sequences.

.. snakemakerule:: cutadapt

Mapping on rRNA
^^^^^^^^^^^^^^^^^^^^^

In order to estimate rRNA rate, a bowtie alignment on ribosomal sequences is performed

.. snakemakerule:: bowtie1_mapping_dynamic

Mapping on reference genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

According to your genome, you can choose to align sequences with bowtie1 or STAR, or both.

Bowtie1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: bowtie1_mapping_dynamic

STAR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: star_mapping

Counting
^^^^^^^^^^^^

FeatureCounts is used for counting

.. snakemakerule:: feature_counts

Reporting
^^^^^^^^^^^^

MultiQC allows to report all bioinformatics tools in a same html file. 

.. snakemakerule:: multiqc
