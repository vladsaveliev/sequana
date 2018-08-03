:Overview: ChIPSeq: Differential expressed peaks analysis
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: BAM, BED and HTML files



Usage
~~~~~~~~~

Example::

    sequana --pipeline chipseq -i data/ -o analysis --no-adapters
    cd analysis
    sbatch snakemake -s chipseq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json --cluster "sbatch --mem={cluster.ram} --cpus-per-task={threads}"

Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/chipseq/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/chipseq/dag.png


Details
~~~~~~~~~

Snakemake ChIP-seq pipelines based on


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/chipseq/config.yaml` to be used with the pipeline.
Each rule used in the pipeline may have a section in the configuration file.
Here are the rules and their developer and user documentation.

