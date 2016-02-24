

Splitting and joining FASTQ files. This could be useful to analyse a large file::

    zcat reads.fastq.gz | split -d -l 2000000 - Block

creates a set of files named Block01, Block02, ... with 0.5M reads (2M lines)
Note that files are not compressed and do not have extension so you will need to fix that::

    for name in Block??
    do
      mv ${name} ${name}.fastq
      gzip ${name}.fastq
    done

You can recombine them::

    zcat Block*.fastq.gz | gzip > allreads.fastq.gw



or to get a randomly selected subset::

    zcat fastq.gz | sort -R | head -10000 > test.fq


split gzipped file::

    gunzip -c bigfile.gz | split -l 400000



