:Overview: Variant calling
:Input: FASTQ file from Illumina Sequencing instrument
:Output: VCF and HTML files
:Config file requirements:
    - samples: file1
    - samples: file2
    - project
    - reference: reference.fasta
Usage
~~~~~~~~~

::

    sequana --pipeline variant_calling --file1 R1.fastq.gz --file2 R2.fastq.gz --project variant --reference reference.fasta
    cd variant
    snakemake -s variant_calling.rules -p --stats stats.txt -j 4

Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/variant_calling/requirements.txt


.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/variant_calling/variant_calling_dag.png

Details
~~~~~~~~

Snakemake variant calling pipeline is based on
`tutorial <https://github.com/ekg/alignment-and-variant-calling-tutorial>`_
written by Erik Garrison. Input reads (paired or single) are mapped using
`bwa <http://bio-bwa.sourceforge.net/>`_ and sorted with
`sambamba-sort <http://lomereiter.github.io/sambamba/docs/sambamba-sort.html>`_.
PCR duplicates are marked with
`sambamba-markdup <http://lomereiter.github.io/sambamba/docs/sambamba-sort.html>`_. 
`Freebayes <https://github.com/ekg/freebayes>`_ is used to detect SNPs and short
INDELs. The INDEL realignment and base quality recalibration are not necessary
with Freebayes. For more information, please refer to a post by Brad Chapman on
`minimal BAM preprocessing methods
<https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/>`_.

The pipeline provides an analysis of the mapping coverage using
`sequana coverage <http://www.biorxiv.org/content/early/2016/12/08/092478>`_.
It detects and characterises automatically low and high genome coverage regions.

Detected variants are annotated with `SnpEff <http://snpeff.sourceforge.net/>`_ if a
GenBank file is provided. The pipeline does the database building automatically.
Although most of the species should be handled automatically, some special cases
such as particular codon table will required edition of the snpeff configuration file.

Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file 
:download:`../sequana/pipelines/variant_calling/config.yaml` to be used with the
pipeline. Each rule used in the pipeline may have a section in the
configuration file. Here are the rules and their developer and user documentation.

Mapping
#########

Add locus in fasta
^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: snpeff_add_locus_in_fasta

BWA
^^^^
.. snakemakerule:: bwa_mem_dynamic

Sambamba markdup
^^^^^^^^^^^^^^^^^^
.. snakemakerule:: sambamba_markdup

Sambamba filter
^^^^^^^^^^^^^^^^^^
.. snakemakerule:: sambamba_filter

Variants Calling
###################

Freebayes
^^^^^^^^^^
.. snakemakerule:: freebayes

Freebayes filter
^^^^^^^^^^^^^^^^^^
.. snakemakerule:: freebayes_vcf_filter

Annotation
####################

.. snakemakerule:: snpeff

Sequana coverage
###################
.. snakemakerule:: sequana_coverage
