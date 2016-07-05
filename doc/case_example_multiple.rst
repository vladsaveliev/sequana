Running multiple pipelines on many files
===========================================

:Description: Running N jobs (N being laqrge) on an LSF cluster using qsub commands


Imagine that you have N pairs of fastq files, where N is large meaning you do
not want to run the pipelines manually. We also want to tell snakemake
to run the the commands on a cluster. Using the **sequana** utility, you can 
type::

    sequana --pipeline quality --glob "*fastq.gz"  
        --cluster "qsub -V -cwd -l mem_total=16G -pe thread 8"

This create N directories with the relevant snakefiles and a script named
**runme.sh**. You will need to run each runme script manually, which is not very
practical so we also add a script named **multirun.sh** to put all required
commands together::

    cd directory1
    sh runme.sh &
    sleep 0.5
    cd ../

    cd directory2
    ...


This should work out of the box but some jobs may fail. If so, it will be
difficult to identify the failed jobs. To help you identifying the jobs, 
we added some utilities in **sequana** in :class:`sequana.misc.FixFailedProjects`.

    from sequana import misc
    fixer = misc.FixFailedProjects()
    fixer()

This will re-generating the a multirun.sh script with only the failed job.



.. note:: For developers, note that although the pipelines are built 
    with Snakemake and that we could have used the **expand** tool to parallelise
    the snakemake across all pairs of files, we we decided not to use that feature. 
    One reason being that this prevents to run a specific rule at will




