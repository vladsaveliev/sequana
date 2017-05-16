.. _installation:

Installation
##########################################


.. _installation_conda:

Installation using Conda
============================


If you have not installed **Sequana**, be aware that it relies on many dependencies
that needs to be compiled (i.e., it is time consumming and requires proper C compilator).
For example, we use Matplotlib, Pandas that requires compilation.
Besides, many pipelines rely on third-party software such as BWA or samtools that are not
Python libraries. Yet, using **conda**, this process is simplied.

Install conda executable
----------------------------

In practice, we do use `Anaconda <https://conda.readthedocs.io/>`_ . We recommend to
install **conda** executable via the manual installer (`download <https//continuum.io/downloads>`). 
You may have the choice
between Python 2 and 3. We recommend to choose a Python version 3.

channels
-----------

When you want to install a new package, you have to use this type of syntax::

    conda install ipython

where **ipython** is the package you wish to install. Note that by default,
**conda** looks on the official Anaconda website (channel). However, there are
many channels available. We will use the **bioconda** channel. To use it, type
these commands (once for all)::

    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda


.. warning:: it is important to add them in this order**, as mentionned on bioconda webpage
(https://bioconda.github.io/).

Create an environment
-------------------------

Once **conda** is installed, open a new shell.
Although this is not required strictly speaking, we would
recomment to create an environment dedicated to Sequana. This environment can
later be removed without affecting your system or conda installation. A
**conda** environment is nothing else than a directory and can be created as
follows::

    conda create --name sequana python=3.5

Then, since you may have several environments, you must activate the **sequana**
environment itself::

    source activate sequana

Here are compulsary packages that must be installed::

    conda install numpy matplotlib pandas snakemake graphviz pygraphviz scipy

Then, depending on the pipelines or standalone applications you want to use,
you will need to install other packages. Here is a list
of dependencies that should be enough to run most of the current
pipelines (commands are split on several lines but you can also
install everything in one go)::

    conda install pysam snpeff biokit bioservices spades khmer pyVCF
    conda install bwa bcftools samtools bedtools picard freebayes fastqc
    conda install kraken krona pigz
    conda install ipython cutadapt jupyter pbr

For atropos::

    pip atropos==1.0.23

or::

    conda install atropos>1.1.2

conda remove tqdm

Install sequana
---------------------

If you have already installed **Sequana** dependencies, this command
should install the latest release posted on http://pypi.python.org/pypi/sequana website::

    pip install sequana --upgrade



.. note:: atropos is an alternative to cutadapt with additional options but same
   type of functionalties and arguments. We use version 1.0.23 and above though. 

.. note:: the denovo_assembly pipelines uses Quast tool, which we ported to
    python 3.5 and was pulled on Quast official github page. This is not
    yet in bioconda but one can it from the quast github (sept 2016). This is
    required for the de-novo pipeline. The denove pipeline also requires GATK, 
    to be installed manually by users (due to licensing restrictions)

.. note:: **Sequana** is not fully compatible with Python 2.7 since a dependency
    (Snakemake) is only available for Python 3.5. However, many core
    functionalities would work under Python 2.7


.. note:: For GATK (variant caller), please go to
   https://software.broadinstitute.org/gatk/download/auth?package=GATK and
   download the file GenomeAnalysisTK-3.7.tar.bz2 ; then type::

    gatk-register GenomeAnalysisTK-3.7.tar.bz2


.. include:: ../docker/README.rst






