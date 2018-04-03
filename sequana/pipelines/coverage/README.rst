:Overview: Coverage
:Input: One BAM or BED file. BED file must have 3 or 4 columns. First column is
    the chromosome/contig name, second column stored positions and third the
    coverage. Fourth optional columns contains a filtered coverage (not used in the
    analysis but shown in the HTML reports)
:Output: multiqc_report.html

Usage
~~~~~~~

The coverage tool takes as input a BAM or a BED file. The BED file must have 3 or 4 columns
as explained in the standalone application (sequana_coverage) `documentation <http://sequana.readthedocs.io/en/master/applications.html?highlight=coverage#sequana-coverage>`_. In short, the first column is the chromosome name, the second column is the position (sorted) and the third column is the coverage (an optional fourth column would contain a coverage signal, which could be high quality coverage for instance).

If you have only BAM files, you can convert them using **bioconvert** tool or
the command::

    samtools depth -aa input.bam > output.bed


If you have a CRAM file::

    samtools view -@ 4 -T reference.fa -b -o out.bam  in.cram

For very large BAM/BED files, we recommend to split the BED file by
chromosomes. For instance for the chromosome  chr1, type::


    # samtools index in.bam
    samtools depth -aa input.bam -r chr1 in.bam > chr1.bed


The standalone or Snakemake application can also take as input your BAM file and
will convert it automatically into a BED file.

Once you have your files in a directory, you can use the standalone
application::

    sequana_coverate --input file.bed


You may also use **Sequanix** in that directory
as follows::

    sequanix -w analysis -i . -p coverage

Go to the second panel, in Input data and then in Input directory. There, you
must modify the pattern (empty field by default meaning search for fastq files)
and set the field to either::

    *.bed

or::

    *.bam


You are ready to go. 

Save the project and press Run. Once done, open the HTML report.


Requirements
~~~~~~~~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/coverage/dag.png


