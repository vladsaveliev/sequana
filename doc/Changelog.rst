Changelog
=============

.. contents::


0.1.15
------------

* CHANGES:

    - coverage: https://github.com/sequana/sequana/issues/302
      add histogram, better stats table. add --output-directory
    - Update docker (add bowtie, subread, firefox)


0.1.14
------------

* CHANGES:

    - fastqc.histogram_sequence_lengths (log2 scale to log10)
    - multi_summary fixed and available for the quality_control pipeline
    - sequana_compressor: add --keep-going option by default so that if a file
      fails, other independent files are processed.
    - snaketools:
          - remove SnakeMakeProfile (not used)
          - remove sequana_check_config (not used)
          - remove deprecated __get_tagname
          - remove ExpandedSnakefile since not required anymore
          - Fix sample_file2 option that was not encoded properly
          - PipelineManager and SequanaConfig use new yaml parser
    - sequana_coverage: -- add back the sample name as prefix of the HTML report
      name -- a BED with two coverage columns is now accepted --
      --download-genbank option added
    - sequana_summary works for the quality_control pipeline
    - Simplify combos of input_directory, input_patter, input_samples, the new
      possible mutually exclusive input parameters of sequana standalone and all
      pipelines.

* BUGS:

    - Kraken: if no reads classified at all, errors were raised and
      quality_control summary report would fail. This is fixed now with a "nodata"
      image being shown.

* NEWS

    - GUI (draft version)
    - fq.gz are now allowed in the pipelines and should be supported in the
      future
    - More tests in particular a ./test/pipelines/ new directory


0.1.13
-----------

* CHANGES:

    - revisited all pipelines so that they can work of multi samples.
    - quality_phix, quqlity and quality_taxon pipelines merged in
      quality_control pipeline
    - running meadian won't fail anymore with odd window size (we add +1)
    - rulegraph is used as well as dag to create figures of the pipelines

* NEWS:

    - compressor: includes dsrc format in addition to bz2 and gz
    - snakemake rule extension for sphinx
    - add a pipeline manager in snaketools to handle all pipelines
    - a designexp module to handle adapter design files


     

0.1.12
-----------

* BUGS:

   - Fix bug in cutadapt pipeline when there is no adapters. Force a dummy
     adapters (XXX) otherwise trimming is performed on read1 only

* NEWS:

    - compressor rule and script available.
    - coverage annotation
    - multiple_summary draft

0.1.11
----------

* NEWS:

   - add a docker
   - sequana_summary standalone
   - sequana_mapping standalone
   - Module has an overview field

* BUG FIXES:

   - cutadapt report handles single-end tables. Fix the reverse complement
     adapter files for the paired-end case

* CHANGES:

    - sequana_standalone: final version with stats



0.1.10 - July 2016
----------------------

* NEWS:

    - sequana_coverage standalone
    - de-novo pipeline

* CHANGES:

    - Remove AdapterDB, a draft version that uses Kraken to detect adapters. Not
      relevant anymore
    - config.yaml is now in each pipeline to have a simplified version
    - sequana can known use single_indexed or multiple_indexed  adapters, which
      are also provided within sequana (Nextera and PCR free cases)
    - Release for production (quality_taxon pipeline)


0.1.7 to 0.1.9 - July 2016
----------------------------

* NEWS:

  - rule data added and used in phix_removal (fastq_sampling + raw data switch)
  - kmer module
  - sequana_taxonomy standalone

* CHANGES:

  - reports are now in ./sequana/reporting
  - MAJOR refactoring of report/ directories in all pipelines to make them
    independent from the temporary analysis, which can then be removed.

* BUGS:

  - Fix running median issue in bedtools (window size larger than contig size)



0.1.6 - June 2016
---------------------


* NEWS:

  - KrakenDownlad class: download kraken_toydv from sequana/data repository or
    minikraken into a local directotry
  - New method in FastQC to show ACGT content
  - Genomecov renamed into GenomeCov
  - Update main script significantly to create multiruns and handle adapters
  - GC content and plot GC vs coverage added in GenomeCov

* CHANGES:

  - sequana_data by default looks into resources/testing directory
  - in fastq module: FastQC a bit faster andFastQRandom class removed
  - add a moving_average function in misc module

* BUGS:

  - sequana_data was showing __init__ and __pycache__ as possible data sets
  - databases: filelist as a list was not implemented
  - in fastq.FastQ extra_head in gzip mode was missing the last row



0.1.5 June 2016
--------------------

* NEWS

  - sequana_taxonomy standalone available (kraken + krona)
  - sequana standalone available
  - quality_taxon pipeline available
  - module coverage for theoretical computations

* CHANGES:

  - module vcf_to_snpeff renamed as snpeff
  - lots of doc added
  - clean adapters module

* BUG:

  - Fix bug in running median (shift)





0.1.4
--------

  - add gallery in the documentation
  - remove old pipelines/ directory
  - sequana standlone refactored (--init option added)
  - Pipeline quality_taxon added
  - Taxonomy rules included
  - Fix the stats image creation


0.1.3
--------

* NEWS
    - Update the pipeline phix_removal


0.1.1 - 0.1.2
--------------
* NEWS
    - Phix pipeline added



0.1.0 April 2006
------------------

* NEWS

  - rules in sequana/rules and pipelines in sequana/pipelines
  - standalone tool called sequana to download a Snakefile and config file
  - modules for post-analysis: bamtools, vcf_filter, fastq, ....
  - Fully tested (90% coverage) and documented
  - Set of reports


