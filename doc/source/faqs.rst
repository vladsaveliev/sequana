
Create a conda environment on IP cluster::

    module load conda
    conda create --name py35 python=3.5
    source condaenvs/py35/bin/activate py35

add channel from where to download packages::

    conda config --add channels r bioconda

Then::

    conda install numpy matplotlib pandas cutadapt pysam bwa bcftools pyvcf samtools snakemake biokit

