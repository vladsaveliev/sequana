:Overview: ATAC-seq: chromatine accessibility pipeline
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: BAM, peaks and HTML files



Usage
~~~~~~~~~

Example::

    sequana --pipeline atac-seq -i data/ -o analysis --no-adapters
    cd analysis
    sbatch snakemake -s atac-seq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json
    --cluster "sbatch --mem={cluster.ram} --cpus-per-task={threads}"

Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/atac-seq/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/atac-seq/dag.png


Details
~~~~~~~~~

Transposons are believed to incorporate preferentially into genomic regions free of nucleosomes (nucleosome-free regions)
or stretches of exposed DNA in general. Thus enrichment of sequences from certain loci in the genome indicates absence
of DNA-binding proteins or nucleosome in the region.
An ATAC-seq experiment will typically produce millions of reads that can be successfully mapped on the reference genome
(with bowtie2). After elimination of duplicates (MarkDuplicates), each sequencing read points to a position on the
genome where one transposition (or cutting) event took place during the experiment. One can then assign a cut count for
each genomic position and create a signal with base-pair resolution (DeepTools Coverage).
Regions of the genome where DNA was accessible during the experiment will contain significantly more sequencing reads
(since that is where the transposase preferentially acts), and form peaks in the ATAC-seq signal that are detectable
with peak calling tools (MACS2).


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


Mapping on the reference genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bowtie2 is used for mapping

.. snakemakerule:: bowtie2_mapping

Identifies duplicate reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating
from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR.

.. snakemakerule:: mark_duplicates


Peak Calling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MACS2 is used for peak calling

.. snakemakerule:: macs2

.. warning:: MACS2 is not included as requirements in Sequana pipelines. We asked for the portage to python 3 (See
https://github.com/taoliu/MACS/issues/179)


DeepTools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DeepTools is a suite of python tools particularly developed for the efficient analysis of high-throughput sequencing
data, such as ChIP-seq, RNA-seq or MNase-seq.

.. snakemakerule:: BamPEfragmentSize
.. snakemakerule:: bamCoverage
.. snakemakerule:: computeMatrix
.. snakemakerule:: multiBamSummary
.. snakemakerule:: plotCorrelation
.. snakemakerule:: plotCoverage
.. snakemakerule:: plotHeatmap


Reporting
^^^^^^^^^^^^

MultiQC allows to report all bioinformatics tools in a same html file.

.. snakemakerule:: multiqc
