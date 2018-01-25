:Overview: Coverage
:Input: One BAM or BED file. BED file must have 3 or 4 columns. First colum is
the chromosome/contig name, second column stored positions and third the
coverage. Fourth optional columns contains a filtered coverage (not used in the
analysis but shown in the HTML reports)
:Output: multiqc_reports.html

Usage
~~~~~~~

First copy the BAM files into a directory. Start **Sequanix** in that directory
as follows (input BAM files must have the extension **.bam**)::

    sequanix -w analysis -i . -p coverage

You are ready to go. 

Save the project and press Run. Once done, open the HTML report.


Requirements
~~~~~~~~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/coverage/dag.png


Details
~~~~~~~~~



Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/coverage/config.yaml` to be used with the pipeline.

