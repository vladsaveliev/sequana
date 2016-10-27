
.. _rules:

Rules
##############

Here below we present a small subset of rules used within the different
pipelines. Rules can be incorporated within a pipeline as follows (given its
name, here e.g. *data*)::

    from sequana import snaketools as sm
    include: sm.module['data']


.. contents::
    :depth: 3

Taxonomy, contaminant
=======================
kraken
------------------
.. snakemakerule:: kraken

fastq_screen
------------------
.. .. snakemakerule:: fastq_screen


fastqc
==========

.. snakemakerule:: fastqc


codecs
===========
The are many codes to be found amongst the rules. For example the rules used by the :ref:`pipeline_compressor` standalone application such as gz_to_fastq:

gz_to_fastq
--------------
.. snakemakerule:: gz_to_fastq

bz2_to_fastq
-------------
.. snakemakerule:: bz2_to_fastq

dsrc_to_fastq
--------------

.. snakemakerule:: dsrc_to_fastq

gz_to_bz2
---------------
.. snakemakerule:: gz_to_bz2

bz2_to_gz
--------------
.. snakemakerule:: bz2_to_gz

bwa_mem
===========
.. snakemakerule:: bwa_mem_dynamic


rulegraph
=============

.. snakemakerule:: rulegraph
.. snakemakerule:: dag


snpeff
===========

.. snakemakerule:: snpeff
.. snakemakerule:: freebayes
.. snakemakerule:: quast

.. sequana/rules/add_locus_in_fasta/add_locus_in_fasta.rules
    sequana/rules/vcf_filter/vcf_filter.rules
    sequana/rules/Misc/conda/conda.rules
    sequana/rules/Misc/cleanup/cleanup.rules
    sequana/rules/samtools_depth/samtools_depth.rules
    sequana/rules/format_contigs/format_contigs.rules
    sequana/rules/digital_normalisation/digital_normalisation.rules
    sequana/rules/fastq_sampling/fastq_sampling.rules
    sequana/rules/samtools_sort/samtools_sort.rules
    sequana/rules/bam_quality_filter/bam_quality_filter.rules
    sequana/rules/syntview/syntview.rules
    sequana/rules/fastq_stats/fastq_stats.rules
    sequana/rules/bwa_index/bwa_index.rules
    sequana/rules/Adapters/skewer/skewer.rules
    sequana/rules/Adapters/alien_trimmer/alien_trimmer.rules
    sequana/rules/Adapters/adapter_removal/adapter_removal.rules
    sequana/rules/Adapters/identify_adapters/identify_adapters.rules
    sequana/rules/Adapters/cutadapt/cutadapt.rules
    sequana/rules/Dev/kraken_adapters/kraken_adapters.rules
    sequana/rules/bwa_bam_to_fastq/bwa_bam_to_fastq.rules
    sequana/rules/indel_realigner/indel_realigner.rules
    sequana/rules/add_read_group/add_read_group.rules
    sequana/rules/create_sequence_dictionary/create_sequence_dictionary.rules
    sequana/rules/mark_duplicates/mark_duplicates.rules
    sequana/rules/spades/spades.rules
