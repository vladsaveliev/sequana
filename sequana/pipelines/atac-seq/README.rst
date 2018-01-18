:Overview: ATAC-seq: chromatine accessibility pipeline
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: BAM, peaks and HTML files



Usage
~~~~~~~~~

Example::

    sequana --pipeline atac-seq -i data/ -o analysis --no-adapters
    cd analysis
    sbatch snakemake -s atac-seq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json --cluster "sbatch --mem={cluster.ram} --cpus-per-task={threads}"

Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/atac-seq/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/atac-seq/dag.svg


Details
~~~~~~~~~

Snakemake atac-seq pipeline s based on a workflow used at Biomics Pole in Institut Pasteur.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/atac-seq/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. Here are the rules and their developer and user documentation.



FastQC
^^^^^^^^^^^

FastQC is used to check quality of sequenced reads.

.. snakemakerule:: fastqc_dynamic


Cutadapt
^^^^^^^^^

Cutadapt is used to trim and filter sequences.

.. snakemakerule:: cutadapt


Bowtie2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bowtie2 is used for mapping

.. snakemakerule:: bowtie2_mapping

MACS2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MACS2 is used for peak calling

.. snakemakerule:: macs2



Reporting
^^^^^^^^^^^^

MultiQC allows to report all bioinformatics tools in a same html file.

.. snakemakerule:: multiqc
