Sequanix Tutorial
====================

Sequanix is a graphical interface dedicated to running Snakemake workflows.

This GUI can be used to load Snakefile and their configuration file. A 
working directory has to be set. Once done, the configuration file can be 
changed in the GUI. Finally, one can run the snakefile and see the progress.
Tooltips are automatically created from the configuration file (if
documented).

Since snakemake has the ability to run jobs locally or on a cluster, this 
application can also be run either locally or a distributed computing
platform (e.g., cluster with slurm scheduler). Of course, this means you can use a X
environment on your cluster (ssh -X should do it)


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
- 

Example: quality_control
------------------------------





























































