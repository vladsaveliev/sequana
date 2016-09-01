
.. _applications:

Applications (standalone)
============================


sequana
---------

:Overview: Creates project(s) to run a **Sequana** pipeline(s)

The **sequana** executable can be used to create pipelines (and associated
config file). For example::

    sequana --pipeline quality --file1 R1.fastq.gz --file2 R2.fastq.gz --project TEST

will create a directory called TEST with a few files such as *quality.rules*,
*config.yaml*, a *runme.sh* and a *README* file.

Valid pipelines can be found using::

    sequana --show-pipelines

There aer many more options and documentation can be obtained using ``--help``
option.

sequana_coverage
--------------------

:Description: Show coverage and interval of confidence to identify under and
    over represtented genomic regions.


sequana_summary
------------------

:Description: Prints basic statistics about a set of NGS input files. Currently
    handles Fastq (gzipped or not) or BED files (coverage).



sequana_taxonomy
--------------------

:Description: Creates a HTML document with Krona and pie chart of taxonomic
    content of a st of FastQ files. Uses Kraken and a dedicated Sequana
    database.

fastq related: fastq_count
-----------------------------

:Description: count number of reads and lines

fastq related: fastq_head
-----------------------------

:Description: Extract head of a fastq files (zipped or not)

