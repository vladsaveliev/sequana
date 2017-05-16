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

There are two kind of dependencies. First, the Python libraries such as
matplotlib or Pandas. Second, the external tools such as BWA (alignment) or
Kraken (taxonomy). The first kind of tools can be installed using Anaconda and the
default conda channel. For instance::

    conda install pandas

The second kind of tools can also be installed using another conda channel
called **bioconda**. For instance::

    conda install bwa

The full list of dependencies will be maintained in the :ref:`installation`
section.


Installation issues
-----------------------


As explained in the previous section, most of the dependencies can be installed
via Conda. If not, pip is recommended. Yet there are still a few dependencies
that needs manual installation. 

quast
~~~~~~~~~

http://quast.bioinf.spbau.ru/manual.html#sec1

::

    wget https://downloads.sourceforge.net/project/quast/quast-4.2.tar.gz
    tar -xzf quast-4.2.tar.gz
    cd quast-4.2

Alternatively, get the source code from their GitHub (takes a while)::

    git clone https://github.com/ablab/quast
    cd quast
    python setup.py install

graphviz
~~~~~~~~~~~~~~~~~~

graphviz provides an executable called **dot**. If you type **dot** in a shell
and get this error message::

    Warning: Could not load
    ...lib/graphviz/libgvplugin_gd.so.6" - file not found

This may be solved by re-installation graphviz using the main anaconda channel
(instead of bioconda)::

    conda install --override-channels -c anaconda graphviz=2.38.0 

:Update April 2017: replace anaconda with conda-forge


matplotlib
~~~~~~~~~~~~~~~~~

If you get errors related to the X connection, you may need to change the
backend of matplotlib. To do so, go in your home directory and in this directory

    /home/user/.config/matplotlib ,

add a file called **matplotlibrc** with the following content::

    backend: Agg

Save, exit the shell, start a new shell.


pysam
~~~~~~~~~~~~~~~~~~

::

    from pysam.libchtslib import *
    ...ImportError: libhts.so.1: cannot open shared object file: No such file or directory


This may be solved by removing conda installation and using pip instead::

     conda remove pysam
     pip install pysam


qt
~~~~~~~~~~~~~~~~~~
::

    from PyQt5.QtWebKitWidgets import QWebView
    ...ImportError: libQt5WebKitWidgets.so.5: cannot open shared object file: No such file or directory

This may be solved by re-installation qt using the main anaconda channel
(instead of bioconda)::

     conda install --override-channels -c anaconda qt


Expected input format
----------------------------

Most of the pipelines and standalone expect FastQ files with the extension
**fastq.gz** meaning that files are gzipped.


Besides, the filename convention is as follows::

    PREFIX_R1_.fastq.gz

that is **_R1_** and **_R2_** indicates the paired or single-ended files and
the PREFIX is used to create directories or reports; it must be present.

.. versionadded:: 0.2
    more flexible tags are now possible in sequana pipelines and sequanix using
    e.g. _R[12] in the **input_readtag** in the configuration file of the
    pipelines.


Sequanix related
----------------------

For question related to Sequanix, we have a dedicated section in
:ref:`sequanix_faqs`.


QXcbConnection issue
----------------------
If you get this error::

    QXcbConnection: Could not connect to display localhost:10.0

this is an issue with your Qt backend. You need to change it to Agg.




















