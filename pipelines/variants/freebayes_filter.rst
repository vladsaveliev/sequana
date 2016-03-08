List of paramater in freebayes VCF
===================================

Primary parameter
------------------

The most important parameter is the QUAL field, which estimates the probability that there is a polymorphism at the loci described by the record.
The value is 1-P(locus is homozygous given the data).

Other parameters
-----------------

- DP: Total read depth at the locus

- RO: Reference allele observation count

- AO: Alternate allele observation count

- SRP: Strand balance Probability for the reference allele

- SAP: Strand balance Probability for the alternate allele

- RPP: Read Placement Probability

- RPPR: Read Placement Probability for reference observation

- EEP: End Placement Probability

- EEPR: End Placement Probability for reference observation


