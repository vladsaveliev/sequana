.. _tutorial:

Tutorial
==========

The following example will show how to run the quality control on a pair of
FastQ files. The data comes from a sequencing (using HiSeq technology) of a
Measles virus. For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.


Quality pipeline
---------------------

**Sequana** comes with standalone applications and pipelines in the form of
Snakefile (`snakemake <https://bitbucket.org/snakemake/snakemake/wiki/Home>`_)

The following example will show how to initialise and run the quality control
pipeline
on a pair of FastQ files.
The data comes from a sequencing (using HiSeq technology) of a
Measles virus. For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.

First, run the sequana standalone application to initialise the pipeline
**quality_control**::

    sequana --pipeline quality_control --working-directory TEST --adapters
    PCRFree

This command downloads the required configuration file(s) in particular
the config file and the pipeline itself. This example should work out of
the box but you may want to look at the
configuration file **config.yaml**. For instance, you may want to change the
reference to the *phix* (by default we use *phix174.fa*, which is provided in
Sequana) or
change the adapter_removal section to your needs (cutadapt parameters, in
particular the forward and reverse complement list of adapters; None by
default).

By default, the output directory is called **analysis** and can be overwritten
with the ``--working-directory`` parameter. Then, run the pipeline and wait for
completion.::

    cd TEST
    snakemake -s quality_control.rules --stats stats.txt -p -j 4 --forceall

The -p option shows the commands, -j 4 means use 4 threads when possible.
Alternatively, there is also a **runme.sh** script.

You should now have a directory with a HTML report corresponding to the sample::

    open Hm2_GTGAAA_L005/report_qc_Hm2_GTGAAA_L005/summary.html


See :ref:`quick_start`


Taxonomy
-------------------------------

Download a toy kraken database designed for this problem (contains only 100
FASTA files mixing measles viruses and others viruses)::


    from sequana import KrakenDownload, sequana_config_path
    kd = KrakenDownload()
    kd.download("toydb")
    database_path = sequana_config_path + "/kraken_toydb"

Then, you may use a Sequana pipeline (see :ref:`pipeline_taxon` and
 :mod:`sequana.kraken`) or this standalone application::

    sequana_taxonomy  --file1 Test_R1.cutadapt.fastq.gz
        --file2 Test_R2.cutadapt.fastq.gz --database  <database_path>

where <database_path> must be replaced with the proper path.


Open the local HTML file krona.html. An example is available
in  `Krona example <_static/krona.html>`_


Variant calling
-------------------

Note that this does the variant calling + snpEff + coverage.
See more information in the :ref:`pipeline_variant_calling` section.



Initialise the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Call **sequana** standalone as follows::

    sequana --pipeline variant_calling --input-directory . --working-directory TUTORIAL

Or use Sequanix. 

Go to the project directory
::

    cd TUTORIAL


Get the genbank reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Assuming the reference is **K01711.1** (Measles virus), we first need to fetch
the genbank file from NCBI::

    from bioservices import EUtils
    eu = EUtils()
    data = eu.EFetch(db="nuccore",id="K01711.1", rettype="gbwithparts", retmode="text")
    with open("measles.gbk", "w") as fout:
        fout.write(data.decode())

Get the FASTA reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will also get the FASTA from ENA::

    from bioservices import ENA
    ena = ENA()
    data = ena.get_data('K01711', 'fasta')
    with open("measles.fa", "w") as fout:
        fout.write(data.decode())


New in v0.10
~~~~~~~~~~~~~~~~

Assuming the genbank and reference have the same name, you can simply
type::

    from sequana.snpeff import download_fasta_and_genbank
    download_fasta_and_genbank("K01711", "measles")

Get a snpEff config file and update it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then you need to initialise a config file for snpEff tool::

    from sequana import snpeff
    v = snpeff.SnpEff("measles.gbk")

Update the snpeff config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit the config file **config.yaml** and add the filename *measles.gbk* in the
snpEff section::

    # snpEff parameter
    snpeff:
        do: yes
        reference: "measles.gbk"

and bwa_ref section::

    # Bwa parameter for reference mapping
    bwa_mem_ref:
      reference: "measles.fa"

.. warning:: In the configuration file, in the mark_duplicates section,
    some output files are huge and requires temporary directory on cluster.


.. warning:: in the configuration file -- coverage section -- note that for short genomes, 
    you may need to decrease the window size.

.. warning:: the mark_duplicates may be changed in the close future to use
   another tool. 


Run the pipeline
~~~~~~~~~~~~~~~~~~~~


::

    snakemake -s variant_calling.rules --stats stats.txt -p -j 4 --forceall


De novo
-------------

The denovo_assembly pipeline can be initialised in the same way::

    sequana --pipeline denovo_assembly --input-directory . --working-directory denovo_test

Go to the **denovo_test** directory and edit the config file. 

.. warning:: this is very time and computationally expensive. The
   **digital_normalisation** section is one that controls the memory footprint.
   In particular, you can check change max-tablesize to a small value for
   test-purposes (set the value to 3e6)




RNA-seq
-------------------


See more information in the :ref:`pipeline_rnaseq` section.
The following example will show how to initialise and run the RNAseq pipeline on a couple of FastQ files.
The data comes from a sequencing (using HiSeq2500 technology) of a saccharomyces cerevisiae strain.
For testing purposes, you can download :download:`Fastq1
<../sequana/resources/data/WT_ATCACG_L001_R1_001.fastq.gz>` and
:download:`Fastq2 <../sequana/resources/data/KO_ATCACG_L001_R1_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.


Initialise the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Call **sequana** standalone as follows::

    sequana --pipeline rnaseq --input-directory . --working-directory EXAMPLE
        --adapter-fwd GATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter-rev GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC

Or use Sequanix.

Go to the project directory
::

    cd EXAMPLE


Get the fasta and GFF reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Assuming the reference is **Saccer3** (Saccharomyces cerevisiae), we first need to fetch
the fasta and the GFF files from SGD to run the pipeline::

    mkdir Saccer3
    cd Saccer3
    wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
    tar -xvzf chromFa.tar.gz
    cat *.fa > Saccer3.fa
    wget http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff -O Saccer3.gff
    rm chr*
    cd ..

.. warning:: All files (fasta, GFF, GTF...) used in RNA-seq pipeline must have the same prefix (Saccer3 in the example)
    and must be placed in a new directory, named as the prefix or not

.. warning:: For the counting step, the RNA-seq pipeline take only GFF files. GTF and SAF files must be integrated soon.

Edit the config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit the config file **config.yaml** and fill the genome section::

    genome:
      do: yes
      genome_directory: ../Saccer3
      name: Saccer3 #path to index name
      fasta_file: ../Saccer3/Saccer3.fa
      gff_file: ../Saccer3/Saccer3.gff
      rRNA_file:
      rRNA_feature: "rRNA"


.. warning:: Note fastq_screen if off by default. It's because sequana not embed a database for this tool.
    If you want to run fastq_screen, please see the manual (https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
    and add the config file in the tool section.


Run the pipeline
~~~~~~~~~~~~~~~~~~~~

On local::

    snakemake -s rnaseq.rules --stats stats.txt -p -j 12 --nolock

on SGE cluster::

    snakemake -s rnaseq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json
    --cluster "qsub -l mem_total={cluster.ram} -pe thread {threads} -cwd -e logs -o logs -V -b y "

on slurm cluster ::

    sbatch snakemake -s rnaseq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json
    --cluster "sbatch --mem={cluster.ram} --cpus-per-task={threads} "

