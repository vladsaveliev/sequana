:Overview: RNASeq
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: 
:Config file requirements:



Usage
~~~~~~~~~

Example::

    sequana --pipeline rnaseq --input-dir .  --output-directory analysis --adapters TruSeq
    cd analysis
    srun  snakemake -s rnaseq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json --cluster "sbatch --mem={cluster.ram} --cpus-per-task={threads}"

Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/rnaseq/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/rnaseq/dag.png


Details
~~~~~~~~~




Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/rnaseq/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. Here are the rules and their developer and user documentation.



FastQC
^^^^^^^^^^^
.. snakemakerule:: fastqc

Fastq_screen
^^^^^^^^^^^
.. snakemakerule:: fastq_screen

Cutadapt
^^^^^^^^^
.. snakemakerule:: cutadapt

Mapping on rRNA
^^^^^^^^^
.. snakemakerule:: bowtie1_mapping_dynamic

Mapping on reference genome
^^^^^^^^^
.. snakemakerule:: bowtie1_mapping_dynamic
.. snakemakerule:: tophat_mapping
.. snakemakerule:: star_mapping

Counting
^^^^^^^^^
.. snakemakerule:: feature_counts

Reporting
^^^^^^^^^^
.. snakemakerule:: multiqc
