:Overview: laa pipeline
:Input: A set of FastQ files or BAM files (see details)
:Output: summary.html

Usage
~~~~~~~

::

    sequana --pipeline laa --input-directory . --working-directory analysis


Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/laa/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/laa/dag.png


Details
~~~~~~~~~

coming soon


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/fastqc/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. 



mutliqc
^^^^^^^^^^^^^^^
.. snakemakerule:: multiqc2

