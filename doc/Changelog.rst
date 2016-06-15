Changelog
=============

.. contents::

0.1.6 June 2016
-----------------

* CHANGES:

  - module vcf_to_snpeff renamed as snpeff
  - lots of doc added
  - clean adapters module

* BUG:

  - Fix bug in running median (shift)


0.1.5 June 2016
--------------------

* NEWS

  - sequana_taxonomy standalone available (kraken + krona)
  - sequana standalone available
  - quality_taxon pipeline available
  - module coverage for theoretical computations
  

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

First functional release.

* NEWS
  - rules in sequana/rules and pipelines in sequana/pipelines 
  - standalone tool called sequana to download a Snakefile + config in a local
directory for further analysis
  - modules for post-analysis: bamtools, vcf_filter, fastq, ....
  - Fully tested (90% coverage) and documented
  - Set of reports

0.0.4 March 2016
-------------------

* NEWS:
    * split samtools module into bamtools and samtools
    * add new class :class:`sequana.bamtools.BAMReport`
* CHANGES
* BUG FIXES

