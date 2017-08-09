:Overview: denovo with canu dedicated to pacbio raw data set followed by quality assessment
:Input: ***.fasta** files  to be found in a directory
:Output: assembly

Usage
~~~~~~~

::

    TODO
    sequana --pipeline pacbio_denovo --input-directory . --working-directory analysis 
    cd analysis
    snakemake -s pacbio_denovo


Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/quality_control/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/pacbio_denovo/dag.png


Details
~~~~~~~~~

- Canu step:
- Busco step:


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/quality_control/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. In the *quality_control* pipeline, we use the *canu*, *busco* rules described here below.


Canu
^^^^^^^^^^^
.. snakemakerule:: canu

Busco
^^^^^^^^^
.. snakemakerule:: busco


