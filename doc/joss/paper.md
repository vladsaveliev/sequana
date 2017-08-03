---
title: 'Sequana': a Set of Snakemake NGS pipelines
tags:
  - snakemake
  - rna-seq
  - variant calling
  - taxonomy
  - denovo
  - pipeline
  - genome coverage  
authors:
 - name: Thomas Cokelaer
   orcid: 0000-0001-6286-1138
   affiliation: 1,2
 - name: Dimitri Desvillechabrol
   orcid: 0000-0002-4947-7876
   affiliation: 2
 - name: Rachel Legendre
   orcid: 
   affiliation: 1,2
 - name: Mélissa Cardon
   orcid:
   affiliation: 2
 affiliations:
 - name: Institut Pasteur - Hub Bioinformatique et Biostatistique - C3BI, USR 3756 IP CNRS
   index: 1
 - name: Institut Pasteur - Biomics Pole
   index: 2
date: 2 August 2017
bibliography: paper.bib
---

# Summary

* Documentation: http://sequana.readthedocs.io
* Repository: https://github.com/sequana/sequana
---
**Sequana** is a Python-based software dedicated to the development of New Generation Sequencing (NGS) pipelines.
The first motivation for this project was to provide NGS pipelines to a sequencing platform (biomics pole at Institut Pasteur, France) that produces more than 200 runs a year combining MiSeq and HiSeq technologies but also long reads technology (PacBio). We decided to use the Snakemake [snakemake:2012] framework to design our pipelines, which ease the decomposition of pipelines into modular sub-units. We currently have 7 pipelines covering quality control, variant calling, long-reads quality, de-novo and RNA-seq analysis (see https://sequana.readthedocs.io for details). Our pipelines are associated with HTML reports based on JINJA templating and Javascript. The reports are used to store the results of a pipelines but also materials required to reprodude the results. **Sequana** is also a Python library that provides tools to perform various analysis tasks (e.g., variant calling filtering). Some of the library components provide original tools that are also available as standalone applications. For instance a fast taxonomic analysis based on Kraken as well as tool to perform exhaustive coverage analysis [coverage:2016].

**Sequana** is an open source project (https://github.com/sequana/sequana). It is developed with the aim
of simplifying the development of new tools (for developers) and the deployment of the pipelines (for users).
The extended documentation (http://sequana.readthedocs.org) and test suite (on Travis.org) provide a high-quality
software that is routinely tested. **Sequana** is now available on bioconda making the installation easier and faster by taking care of the dependencies (e.g., samtools, bwa standalones, or Python librairies). 

Finally, for end-users, we also developed a Graphical interface called **Sequanix** [sequanix:2017] developed with the PyQt framework. **Sequanix** standalone exposes all **Sequana** pipelines (Snakemake pipelines) within an easy-to-use interface. Within the graphical interface, the configuration file used by Snakemake are automatically loaded and can be edited by end-users with dedicated widgets. We made the interface generic enough that not only Sequana pipelines can be run interactively but also any Snakemake pipelines.


![](sequana.png)

# Future works

**Sequana** is an on-going project. Although the project has reached a mature stage with stable pipelines, new pipelines will be including on demand or based on new technologies. 

# References
