Example of a snakemake dedicated to variant calling.
=====================================================

This pipeline does the variant calling of a bam file with freebayes.
Before the calling, GATK best practices are done and mark duplicate was optionnal.
The vcf from freebayes is filtered by a python module.

TODO:
-------
- SNPeff
