.. _installation:

Installation
##########################################

In order to use Sequana, you can either install from source (or conda), or use one of the 
available container (Docker or Singularity).


Quick installation
=====================

* Sequana is available on conda/bioconda as pre-compiled package (**conda install sequana**).
* If you prefer to install everything yourself, the source code is available on
  github (http://github.com/sequana/sequana) and releases of the source code are posted on Pypi (**pip install sequana**). 


.. _installation_conda:


From bioconda (Recommended !! )
-------------------------------------

Sequana is on `bioconda <https://bioconda.github.io/>`_. You can follow these `instructions <http://bioconda.github.io/recipes/sequana/README.html>`_ or type::

    conda install sequana

You would need these special channels::

    conda config --add channels r
    conda config --add channels default
    conda config --add channels conda-forge
    conda config --add channels bioconda


If you have already set the channels, please check that the order is correct.
With the following command::

    conda config --get channels

You should se::

    --add channels 'r'   # lowest priority
    --add channels 'defaults'
    --add channels 'conda-forge'
    --add channels 'bioconda'   # highest priority



This does not provide all dependencies needed by the different pipelines. So,
you may need to install extra packages as listed in this requirement file that
can be used with conda::

    conda install --file https://raw.githubusercontent.com/sequana/sequana/master/requirements_pipelines.txt





From Pypi website (released source code)
------------------------------------------
If you do not want to use **conda**, we provide releases on the Python Package Index website (pip tool)::

    pip install sequana
    pip install PyQt5

See below for dependencies.

From GitHub Source code
------------------------------

Finally, if you are a developer, you can install **sequana** from source::

    git clone git@github.com:sequana/sequana.git
    cd sequana
    python setup.py install

See below for dependencies.

Dependencies
=================

With the methods above, you should have a working libary of Sequana and its
standalones (e.g. Sequanix). However, if you wish to use all the pipelines,
additional tools and libraries are required. Some are available on PyPi (Python
software) but others will be available only on BioConda.


Installation using Conda
============================

If you have not installed **Sequana**, be aware that it relies on many dependencies
that needs to be compiled (i.e., it is time consumming and requires proper C compilator).
For example, we use Matplotlib that requires compilation.
Besides, many pipelines rely on third-party software such as BWA or samtools that are not
Python libraries. Yet, using **conda**, this process is simplified.

Install conda executable
----------------------------

In practice, we do use `Anaconda <https://conda.readthedocs.io/>`_ . We recommend to
install **conda** executable via the manual installer (`download <https//continuum.io/downloads>`_). 
You may have the choice
between Python 2 and 3. We recommend to choose a Python version 3.

Add conda channels
---------------------

When you want to install a new package, you have to use this type of syntax::

    conda install ipython

where **ipython** is the package you wish to install. Note that by default,
**conda** looks on the official Anaconda website (channel). However, there are
many channels available. We will use the **bioconda** channel. To use it, type
these commands (once for all)::

    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

.. warning:: **it is important to add them in this order**, as mentionned on bioconda webpage
    (https://bioconda.github.io/).

Create an environment
-------------------------

Once **conda** is installed, open a new shell.
Although this is not required strictly speaking, we would
recommend to create an environment dedicated to Sequana. This environment can
later be removed without affecting your system or conda installation. A
**conda** environment is nothing else than a directory and can be created as
follows::

    conda create --name sequana_env python=3.5

Then, since you may have several environments, you must activate the **sequana**
environment itself::

    source activate sequana_env

Install sequana via conda (bioconda)
-------------------------------------

Finally, just type::

    conda install sequana

This should install most of the required dependencies. However, you may need to
install more packages depending on the pipeline used. To install all required
packages, you may use this command::

    conda install --file https://raw.githubusercontent.com/sequana/sequana/master/requirements.txt
    conda install --file https://raw.githubusercontent.com/sequana/sequana/master/requirements_pipelines.txt

Developers, can also install pytest and sphinx::

    conda install --file https://raw.githubusercontent.com/sequana/sequana/master/requirements_dev.txt

We would also recommend those tools::

    conda install ipython jupyter 

.. note:: atropos is an alternative to cutadapt with additional options but same
   type of functionalties and arguments. We use version 1.0.23 and above though. 

.. note:: the denovo_assembly pipelines uses Quast tool, which we ported to
    python 3.5 and was pulled on Quast official github page. This is not
    yet in bioconda but one can get it from the quast github (sept 2016). This is
    required for the de-novo pipeline. The denovo pipeline also requires GATK, 
    to be installed manually by users (due to licensing restrictions)

.. note:: For GATK (variant caller), please go to
   https://software.broadinstitute.org/gatk/download/auth?package=GATK and
   download the file GenomeAnalysisTK-3.7.tar.bz2 ; then type::

    gatk-register GenomeAnalysisTK-3.7.tar.bz2


singularity
================

Install singularity (http://singularity.lbl.gov/). Download our image::

    singularity run shub://sequana/sequana:master

and use it. For instance, to use sequana_coverage executable::

    singularity exec sequana-sequana-master.img sequana_coverage --help

or sequanix::

    singularity exec sequana-sequana-master.img sequanix



.. include:: ../docker/README.rst






