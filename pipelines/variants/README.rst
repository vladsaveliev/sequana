Example of a snakemake dedicated to variant calling.
=====================================================

This pipeline does the variant calling of a bam file. The first rule is markDuplicate. It is an optional rule. You must indicate in the config file if you wish to activate this rule. (true/false)

The second rule indexes your bam file. And finally, the last rule call freebayes to detect snp and indel variants.

TODO:
-------
- Mapping rule ?

- Refine freebayes params

- VCF filter.

- SNPeff
