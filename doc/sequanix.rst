Sequanix Tutorial
====================

.. contents::
   :local:

Introduction
---------------
**Sequanix** is a graphical user interface (GUI) that can be used to run Snakemake workflows.
The standalone name is **sequanix** (small caps) and is part of **Sequana**
library.

The primary goal was to provide a GUI to easily run **Sequana pipelines** (designed
as Snakemake workflows).

However, we extended the interface so that it can handle other Snakemake
workflows, referred to as **Generic pipelines** in the GUI.


.. figure:: _static/sequanix/sequanix.png
    :width: 80%

    Snapshot of the Sequanix graphical user interface (GUI)

In this section, we first show how to run one of our Sequana pipeline (quality
control pipeline). Second, we show how to run a **Generic pipeline** that is not
part of Sequana. For these two examples, the computation is done locally.
However, one strength of Snakemake pipelines is that they can be executed on
various cluster without changing the pipeline itself. The
third section shows how to run the analysis on a cluster (SLURM and SGE frameworks).


Snakemake pipelines are made of 2 parts: a pipeline and an optional
configuration file; The pipeline may be called **Snakefile**. It contains the code
of the pipeline itself. Keep in mind that
in the Snakefile, developer may link the pipeline to an external configuration
file: the **config** file, which is encoded in :term:`YAML` or JSON format.


Sequana pipeline: the quality control example
----------------------------------------------------

For all Sequana pipelines, we need:

#. To select the pipeline by its name (the configuration file is loaded automatically)
#. A directory where to find the data
#. A working directory where analysis will be run and results stored

Prerequisites: get some data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following example will show how to run a quality control pipeline on a pair of FastQ files. The data comes from a sequencing platform (using HiSeq technology) of a Measles virus. For testing purposes, you can download :download:`R1 <../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and :download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.

Save the file in a directory (let us call it testing) and move to that directory
(e.g. *cd testing*). Starts **sequanix** in that directory.



Select the quality control pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First you need to select the pipeline of interest (here the quality_control).
In the following figure, you  need to select the sequana pipelines tab (arrow
1), then the pipeline section (arrow 2) and finally the pipeline in the dropdown
box (arrow 3)

.. figure:: _static/sequanix/sequanix_tuto_qc_selection.png

Once done, the configuration file of the pipeline will be loaded in the **Config
parameters** tab (arrow 4).

Fine tune the config parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning:: each Sequana pipeline requires expertise. Check the online
   documentation for help (:ref:`tutorial`, :ref:`pipelines`).

One major interest of **Sequanix** is that the Snakemake configuration file is
loaded and can then be changed dynamically. The advantages of the config
parameters section are diverse:

- one does not need to use an esoteric editor to change the configuration
- file can be selecting thanks to a file browser. If no file is selected, the button is red (green otherwise)
- Some button have dedicated widgets (e.g. in the figure above, the number of threads has its own dropbox limited typing errors)
- Boolean have their own checked button
- etc


.. note:: For developers, please see the :ref:`developer` section.

Save, check and run the project
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the parameters have been selected, one can save the project (arrow 5)







This GUI can be used to load Snakefile and their configuration file. A
working directory has to be set. Once done, the configuration file can be
changed in the GUI. Finally, one can run the snakefile and see the progress.
Tooltips are automatically created from the configuration file (if
documented).

Since snakemake has the ability to run jobs locally or on a cluster, this
application can also be run either locally or a distributed computing
platform (e.g., cluster with slurm scheduler). Of course, this means you can use a X
environment on your cluster (ssh -X should do it)



Preference dialog
---------------------

.. figure:: _static/sequanix/preferences_dialog.png
   :scale: 80%

   Preferences dialog. This dialog is accessible via the menu or the short Ctrl+P.
   It contains general options to tune Sequanix's behaviour


Snakemake dialog
--------------------

The Snakemake dialog contains 3 sub tab: the local, cluster and general tabs.

.. figure:: _static/sequanix/snakemake_dialog_local.png
    :scale: 80%

    The **local** tab contains only one option to set the number of local cores
    to be used. By default it is the number of available cores on the machine
    used.

.. figure:: _static/sequanix/snakemake_dialog.png
    :scale: 80%

    The tab **cluster** contains parameters related to the execution of the Snakemake pipeline can be set (e.g. specific job scheduler information or number of CPUs to be used). 


.. figure:: _static/sequanix/snakemake_general.png
    :scale: 80%

    In the General tab, check boxes related to Snakemake are available. Any other options can be set in the editable line at the bottom. 


config example
~~~~~~~~~~~~~~~~~~~
Snakemake pipeline are usually associated with a configuration file. For instance all Sequana pipelines have their own configuration file named \textit{config.yaml}.  Although Snakemake allows the configuration to be in YAML or JSON format, we recommend to use YAML only. An example is shown in \ref{lst:listing}. In the example, the YAML is made of only one section ( bwa\_mem ). That section contains 4 arguments. This type of file can be read by Snakemake pipeline and sections are exposed as a Python dictionary in the pipeline namespace. YAML files are human-readable and can also be commented. Comments can be added anywhere but must start with the \# sign.

In Sequanix, comments after section or arguments are ignored. However, we encourage developers to have a self-content documentation before the section. The comments should be encoded with the reStructuredText syntax (http://docutils.sourceforge.net/rst.html). It is an easy-to-read plaintext markup syntax and parser system. It is used in Python library to document codes but also in various software projects, generally via Sphinx syntax, which is an augmented version of reStructuredText. If comments are properly encoded, then the text is extracted and interpreted in Sequana. Finally, a tooltip showing the corresponding HTML code is shown in the Sequanix interface as soon as the user moves the mouse cursor on a section. The Fig.~\ref{fig:tooltip} shows how the listing \ref{lst:listing} is transformed and shown as a tool tip. 


First, run the sequana standalone application to initialise the pipeline **quality_control**::



config example
-----------------

::

    ############################################################
    # BWA used to remove a contaminant
    #
    # :Parameters:
    #
    # - do: if unchecked, this rule is ignored
    # - reference_file: the name of the reference file to be found
    #        in the analysis directory.
    # - index_algorithm: the BWA index algorithm
    # - options: any options recognised by BWA tool
    # - threads: number of threads to be used
    #
    bwa_mem:
        reference_file: "phiX174.fa"
        index_algorithm: "is"
        options: '-T 30'
        threads: 4



prerequisites
------------------

Sequanix allows users to select Sequana pipeline, set configuration files
interactively and run the pipeline using Snakemake behind the scene. 

The motivation is to expose complex pipelines via a simple graphical interface.

So, before using Sequanix you must know what pipeline you want to use.

Sequanix
-----------

installation: cf installation de sequana

taper sequanix

!! ou demarrer sequanix ? wherever but we recommend to start it where the data
is.

:Cluster usage: -X 


Running analysis
-------------------

Snapshot sequanix avec les sections mises en valeurs. Entour en rouge les
sections avec labels I, II, III, IV. C'est la figure 1 de l'article.

pipelines are defined by
- a Snakemake file, which describes the pipeline itself
- a config file, where users can fine-tune the pipeline (this may be optional)
- a working directory where we save the pipeline / config and run the analysis
(in general)
- information about input data set (files, directory) 

The first step is to define a project. This is done in window I

Step 1
~~~~~~~~
I. Select a pipeline
I.A Sequana pipeline
...


Step 2: configuration via the form
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Snapshot 

- ability to switch on/off some rules/tools/steps 
- dropdown widget
- file browser to be filled

we cannot details the config for the pipelines so we need to refer to the
pipeline section again.

Once happy we the config, SAVE the project. This copies the config and pipeline
files into the WORKING DIRECTORY. If exists already what's going on ?

Step 3: run
~~~~~~~~~~~~~~~~

- Run. describe progress bar 
- Stop
- Unlock
- Save --> enable the RUN


Section II
~~~~~~~~~~~~~~~
- local/cluster
- if cluster --> preferences 

With Slurm::

    srun --x11 sequanix

FAQS
---------

- What to do if the RUN fails
- 

Browser
----------
- motivation and limitations

others
-----------
- import config local
- start sequanix with options::

   sequanix -i . -p quality_control -w analysis

- Generic pipeline can re-use widgets using _file and other semantic.


Generic example
------------------

For any other Snakemake workflows, we need:

#. To select a Snakefile (extension .rules)
#. To select a configuration file (optional)
#. A working directory where analysis will be run and results stored

























































