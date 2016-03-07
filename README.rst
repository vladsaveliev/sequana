SEQUANA
############


>>>>>>>>>> WIP <<<<<<<<<<<<<<


**Sequana** includes a set of pipelines related to NGS (new generation sequencing). 

It will provide a set of modular pipelines and reports associated to them.


Installation
=================


::

    pip install sequana


Some dependencies required include matplotlib, pandas, cutadapt, pysam. If you
are new or starting with Python, we strongly recommand to use anaconda and to
install those dependencies::

    conda install matplotlib pandas cutadapt pysam

although the code is Python2.7 and Python3.5 compatible, a dependency
(Snakemake) only supports Python3.5 for the moment so, we will assume you have a
Python3.5 version installed.


Example
==========

::

    # contaminant and adapter removal search
    from sequana import biomics


    pipeline = biomics.Biomics("*fastq.gz") # assumm two paired-end fastq files
    pipeline.create()
    pipeline.start()
    pipeline.report()

Behind the scene, a Snakefile and a config file are created in the local
directory. The snakefile can be executed indepedently::

    snakemake Snakefile -p




