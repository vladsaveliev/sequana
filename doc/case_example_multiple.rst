Running multiple pipelines on many files
===========================================


Imagine that you have N pairs of fastq files, where N is large meaning you do
not want to run the pipelines manually. With the **sequana** utility, you can 
first run::

    sequana --pipeline quality --glob "*fastq.gz" 

This create N directories with the relevant snakefiles and a script named
**runme.sh**

In addition, a multirun.sh is created


Jobs may fail. In which case, you will need to figure out the directories where
the scripts faile. This is just not feasible so, we also provide an utility to
create the multirun file with only the jobs that failed. ::

    from sequana import misc
    fixer = misc.FixFailedProjects()
    fixer()


.. note:: For developers, note that although the pipelines are built 
    with Snakemake and that we could have used the **expand** tool to parallelise
    the snakemake across all pairs of files, we we decided not to use that feature. 
    One reason being that this prevents to run a specific rule at will




