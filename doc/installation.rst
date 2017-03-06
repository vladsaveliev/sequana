Installation
##########################################


.. _installation_conda:
Installation using Conda
============================

If you have already installed **Sequana** dependencies, this command
should install the latest release posted on http://pypi.python.org/pypi/sequana website::

    pip install sequana --upgrade

If you have not install **Sequana**, be aware that it relies on many dependencies
that needs to be compiled (i.e., it is time consumming and requires proper C compilator).
For example, we use Matplotlib, Pandas that requires compilation.
Besides, many pipelines rely on third-party software such as BWA or samtools that are not
Python libraries. 

In practice, we do use `Anaconda <https://anaconda.org/>`_ and in
particular the **bioconda** channel, which can be
added to your environment as follows (once Anaconda is installed)::

    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda

**It is important to add them in this order**. 


This is the recommended way and we will not support manual installation.

Here are compulsary packages that must be installed::

    conda install numpy matplotlib pandas snakemake graphviz scipy

Then, depending on the pipelines or standalone applications you want to use,
you will need to install other packages. Here is a list
of dependencies that should be enough to run most of the current
pipelines (commands are split on several lines but you can also
install everything in one go)::

    conda install pysam snpeff biokit bioservices spades khmer pyVCF
    conda install bwa bcftools samtools bedtools picard freebayes fastqc
    conda install kraken krona pigz


.. note:: the denovo_assembly pipelines uses Quast tool, which we ported to
    python 3.5 and was pulled on Quast official github page. This is not
    yet in bioconda but one can it from the quast github (sept 2016). This is
    required for the de-novo pipeline

.. note:: **Sequana** is not fully compatible with Python 2.7 since a dependency
    (Snakemake) is only available for Python 3.5. However, many core
    functionalities would work under Python 2.7


.. include:: ../docker/README.rst



