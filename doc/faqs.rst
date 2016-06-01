FAQS
======

Conda related
---------------

Create a conda environment on IP cluster::

    module load conda
    conda create --name py35 python=3.5
    source condaenvs/py35/bin/activate py35

add channel from where to download packages::

    conda config --add channels r bioconda


What are the dependencies
-----------------------------
There are quite a few dependencies and depends on the pipelines. This no
exhaustive list contains all the python libraries used by default::

    conda install numpy matplotlib pandas pysam bwa samtools snakemake biokit

However, other packages may be required ::

    conda install cutadapt bcftools pyvcf
