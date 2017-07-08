Changelog
=============

.. contents::

2017
------

0.3 April-June 2017
~~~~~~~~~~~~~~~~~~~~~~~~~~~~é

* BUG FIXES:

    - sequanix:
        - rulegraph issue on SLURM system. Avoid the os.chdir
    - fastq_samples/ fastq module: fix histogram_gc_content maximum range
    - rulegraph rule: fix issue #405 (spaces in path to snakefile)
    - genome coverage was buggy for multi chromosome and circular option on. Fixed
    - adapters/expdesign modules: fixe the case of design files with same sample
      name and same index but different lanes.
    - sequana_coverage. Fix Issue #416 (float division by zero)

* CHANGES:

    - sequanix:
        - snakemake output is now cleared
    - quality_control pipeline: default to atropos instead of cutadapt for
          adapter trimming. Kraken: remove classified reads and keep
          unclassified. Unclassified reads are now compressed.
          unclassified reads that are also compressed now.

* NEW:

    - pacbio module: cleanup and add funcion to convert input BAM into Fasta
    - sequence module: Repeats class added
    - new Snakemake pipeline called qc_pacbio to perform quick QC and taxonomy analysis
          for pacbio
    - add ORD, CDS, GC SKEW in sequence module.


0.2. - March - April 2017
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* NEWS:
    - RNA-seq pipeline added (single-end only, paired-end upcoming)
      including all indexes for RNA-seq
    - Hierarchical kraken available
    - add new standalone called **sequana_fox** to expose the pyqt5 browser.
    - Sequanix first release
    - final version of the variant calling, denovo, quality_control and rna-seq
      pipelines.

* CHANGES:

    - Sequanix/Sequana:
      - config file can have the yml extension (in addition to yaml)
      - dropdown widgets in the form based on the docstrings in the config file
      - can import config to override default sequana config file
      - subprocesses killed when the main pipeline is stopped

0.1.21 - Feb 2017
~~~~~~~~~~~~~~~~~~~~~~~~

* NEWS:

    - add sequana_debug_level function at top level to switch verbosity of
      informative messages (default is WARNING).
    - add pacbio module  #351
    - quality control pipeline: atropos can be used in place of cutadapt #346

* CHANGES:

   - Running Median is 10 times faster #345
   - sequana_coverage:  (1) --file1 alone was not working (2) automatically copy
     cluster-config in working directory and update runme.sh accordingly #342
   - sequana standalone:
       - handles cluster_config Snakemake option
       - add error message when adapter name is incorrect
   - sequanix: the help dialog is now created inside designer and has a proper
     scrollable browser dialog. cluster_config Snakemake option is also handle.
   - Remove galleria JS lib and related files (htmltools)
   - sequana_coverage: add --logging-level option

* BUG:

    - Fix #352 : allow gc window size to be even (warning is shown and +1 to
      window size)
    - Fix # 354: cutadapt report that was mixing up R1/R2 trimming in the images.
    - --output-directory in sequana_coverage was failing 
    - in coverage, centralness was buggy (regression) and use number of ROIs
      instead of the total base length #347
    - Fix multi_report summary for single end case #349

0.1.20 - Feb 2017
~~~~~~~~~~~~~~~~~~~~~~~~

* CHANGES: 

    - remove pyquickhelper dependencies and add a simple rest2html function in
      misc module.

0.1.19 - Feb 2017
~~~~~~~~~~~~~~~~~~~~~~~~


* CHANGES:

    - misc module: factorise on_cluster() function used in compressor scripts to
        be used in other tools such as sequanix
    - compressor: limits max number of jobs to 20 (can be bypass manually),
      prevent run on TARS if snakemake-cluster not provided. 
    - rules:
        - dag: now the snakemake is called inside a temporary directory to avoid
          clash with the current snakemake process. This avoid error message. 
          Fixes https://github.com/sequana/sequana/issues/331
    - __init__ was optimized as well as many modules to make use of the lazy
      import mechanism. The reporting package is not part of the exposed module. 
      So::

         from sequana import BAMReport

      is now::

         from sequana.reporting.report_bam import BAMReport

* NEWS:

    - Sequanix stable version
    - add TrueSeq adaptors
    - add lazy import mechanism to speed up the time to import sequana, which 
      speeds up the --help in the standalone 


0.1.17/0.1.18 - Jan 2017
~~~~~~~~~~~~~~~~~~~~~~~~

:Main NEWS: The GUI was completed and the current pipelines stabilised (RNA-seq,
    quality control, variant calling). The test suite was switched from nosetests to
    pytest, in particular to perform tests more eaasily on the Qt GUI. 


* BUG Fixes:

    - experimental design and adapters API simplified fixing a few bugs in the
      process. Doc and tested finalised.
    - Fix cutadapt rules, which was not filling the fwd and rev properly anymore
      when using the design file.
    - in sequana main script, --reference was used by quality_pipeline only.
      Now, available for all.
    - Fix the main script for the reference in variant calling pipeline.


* CHANGES:

    - sequana_compressor: for conversion from e.g gz to bz2, use a pipe instead
        of double IO. Updated docs and tests ready for production.
    - sequana standalone: 
      - --pattern changed to --input-pattern
      - --output-directory changed to --working-directory
    - remove pipetools module (obsolet)
    - GUI revisited with qt designer + can now also read any snakefile/config
      file combo (not just sequana pipelines)
    - RULES: adapters can now use adapter_type without a design (fwd and rev
      gets filled automatically)

* NEWS:

    - add rubicon adapters
    - add ability to read JSON in SequanaConfig

2016
----------
0.1.16
~~~~~~~~~~~

* BUG Fixes:

    - Fix sequana_taxonomy (https://github.com/sequana/sequana/issues/308)
    - Fix typo in sequana_coverage for multiple chromosome (https://github.com/sequana/sequana/issues/307)

* NEWs:

    - SequanaConfig can read back a SequanaConfig instance
    - Added a DummyManager for minimalist manager to create reports


0.1.15
~~~~~~~~~~~

* CHANGES:

    - coverage: https://github.com/sequana/sequana/issues/302
      add histogram, better stats table. add --output-directory
    - Update docker (add bowtie, subread, firefox)
    - snaketools:
          - empty strings are kept as empty strings (not None)
          - remove check() method in SequanaConfig
          - cleanup (removing of templates) ca be switch off

0.1.14
~~~~~~~~~~~

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
~~~~~~~~~~~

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
~~~~~~~~~~~

* BUGS:

   - Fix bug in cutadapt pipeline when there is no adapters. Force a dummy
     adapters (XXX) otherwise trimming is performed on read1 only

* NEWS:

    - compressor rule and script available.
    - coverage annotation
    - multiple_summary draft

0.1.11
~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~


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



prior 0.1.5 June 2016
~~~~~~~~~~~~~~~~~~~~~~

* NEWS

  - sequana_taxonomy standalone available (kraken + krona)
  - sequana standalone available
  - quality_taxon pipeline available
  - module coverage for theoretical computations
  - add gallery in the documentation

* CHANGES:

  - module vcf_to_snpeff renamed as snpeff

* BUG:

  - Fix bug in running median (shift)

