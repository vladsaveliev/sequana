.. _tutorial:

Tutorial
==========

Following the introductory example in :ref:`quick_start`, we will look at other pipelines such as
the taxonomic classification, variant calling and coverage. For completeness,
reproduce the quality pipeline here below. 


The following example will show how to run the quality control on a pair of
FastQ files. The data comes from a sequencing (using HiSeq technology) of a
Measles virus. For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory. 


Quality pipeline
---------------------

.. todo:: WIP


Taxonomy
-------------------------------

Download a toy kraken database designed for this problem (contains only 100
FASTA files mixing measles viruses and others viruses)::


    from sequana import KrakenDownload, sequana_config_path
    kd = KrakenDownload()
    kd.download("toydb")
    database_path = sequana_config_path + "/kraken_toydb"

Then, you may use a Sequana pipeline (see :mod:`sequana.kraken`) or this standalone 
application::

    sequana_taxonomy  --file1 Test_R1.cutadapt.fastq.gz --file2 Test_R2.cutadapt.fastq.gz --database
        <database_path>

where <database_path> must be replaced with the proper path.


Open the local HTML file krona.html. An example is available 
in  `Krona example <_static/krona.html>`_


Variant calling
-------------------

Note that this does the variant calling + snpEff + coverage. 
See more information in the :ref:`pipeline_variant_calling` section.



Initialise the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~




:: 

    sequana --init variant_calling --input-dir . --project TUTORIAL
    cp variant_calling.rules Snakefile # this will not be required in the future




Get the genbank reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming the reference is **K01711.1** (Measles virus), we first need to fetch
the genbank file rfom NCBI::

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

Get a snpEff config file and update it 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then you need to initialise a config file for snpEff tool::

    from sequana import snpeff
    v = snpeff.SnpEff("measles.gbk")


Edit the config file **config.yaml** and add the filename *measles.gbk* in the
snpEff section::

    # snpEff parameter
    snpeff:
        do: yes
        reference: "measles.gbk" 

and bwa_ref section:: 

    # Bwa parameter for reference mapping
    bwa_ref:
      reference: "measles.fa"



Run the pipeline
~~~~~~~~~~~~~~~~~~~~


::

    snakemake -s Snakefile --stats stats.txt -p -j 4 --forceall
