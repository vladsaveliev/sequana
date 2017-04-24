Overview
############

.. contents::

**Sequana** provides standalone applications (e.g., **sequana_coverage**,
**sequana_taxonomy**) and pipelines in the form of Snakefiles. Although the standalone
applications are usually simpler, they may not have all features or parameters
offered by the pipelines.

The :ref:`Tutorial`, :ref:`pipelines`, :ref:`Gallery` and :ref:`case_examples` 
sections provide many examples on their usage. 


This section will not describe all available standalones and pipelines.
We will focus on one example (coverage) to show how one can use
the **Sequana** library, or standalone application, or pipeline to get
information about the coverage of a set of mapped reads onto a reference.


**Sequana** library
========================

Example 1 : running median on coverage
----------------------------------------

**Sequana** is a Python library. It contains many functionalities, which are
fully documented and available in the :ref:`references` section. We can first
look at the coverage contained within a BED file using the library. First, we
need some data. **Sequana** provides some test examples, which can be accessed
using :func:`~sequana.sequana_data` function. The test case is a virus (about
18,000 bases)::

    from sequana import sequana_data
    filename = sequana_data('JB409847.bed')


We can then use the :class:`~sequana.bedtools.GenomeCov` class to read the
file::

    from sequana import GenomeCov
    gc = GenomeCov(filename)

Select a chromosome (first one) and compute the running median::

    chrom = gc[0] 
    chrom.running_median(n=5001, circular=True)
    chrom.compute_zscore()

and finally plot the coverage together with confidence interval (3 sigma)::

    chrom.plot_coverage()


.. plot::

    from sequana import sequana_data
    filename = sequana_data('JB409847.bed')
    from sequana import GenomeCov
    gc = GenomeCov(filename)

    chrom = gc[0]
    chrom.running_median(n=4001, circular=True)
    chrom.compute_zscore()
    chrom.plot_coverage()


Example2: read a fastq file
------------------------------

Let us use the :class:`FastQC` class to get the distribution of the bases ACGT
across all reads of a FastQ file.


.. plot::

    from sequana import FastQC
    from sequana import sequana_data
    filename = sequana_data("test.fastq")

    fastqc = FastQC(filename)
    print(fastqc.fastq)
    for x in 'ACGT': 
        fastqc.get_actg_content()[x].hist(alpha=0.5, label=x, histtype='step', lw=3, bins=10)

    from pylab import legend
    legend()



Many more functionalities are available. The reference guide should help you.

**Sequana** standalones
=========================

The Python example about the coverage is actually quite useful. We 
therefore decided to provide a standalone
application. There are other standalone applications listed in
:ref:`applications` section.

The one related to the coverage example shown above is named
**sequana_coverage**. If you have a BED file, type::

    sequana_coverage  -i <BEDFILENAME> 

If your organism has a circular DNA, add ``-o``. You can play with the window
size for the running median using ``-w``.

Using the BED file and reference mentionned in the previous section you should
obtain the same figure as above.

An additional feature is the report using  ``--show-html`` option.

**Sequana** pipelines
=======================

In **Sequana**, in addition to the library and standalone applications, we also
provide a set of pipelines (see :ref:`pipelines` section). The coverage tools
described so far do not have a dedicated pipeline but is part of a more general
pipeline called :ref:`pipeline_variant_calling`. Instead of describing in
details that pipeline, let us explain the way pipelines can be created and run.

Manually
------------

Pipelines are made of a Snakefile (a Makefile using Python) and an associated
config file. Pipelines can be downloaded from the **Sequana** 
`pipeline directory <https://github.com/sequana/sequana/tree/master/sequana/pipelines>`_
as well as the config file named **config.yaml**.

Copy the pipeline (ending in .rules) and the configuration file in a local
directory. The config file is a generic template file and some fields must be
changed. For instance the beginning of the file looks like::

    # list of your input file
    samples:
        file1: "%(file1)s"
        file2: "%(file2)s"

For pipelines that takes FastQ files as inputs, the string **%(file1)s** must be 
replaced by a valid filename. If you do not have a second file, remove the next
line (file2). Other similar fields must be filled if required by the pipeline.

Then, a pipeline must be executed using the executable **snakemake**. If you
choose the **variant_calling** pipeline, the file is executed as follows::

    snakemake -s variant_calling.rules

This will search for the **config.yaml** file locally. One good feature is that
if you interrupt the pipeline (or if it fails), you can fix the problem and
re-run the command above without executing the parts of the pipelines that were
succesfully run. If you want to start from scratch, add ``--forceall`` option::

    snakemake -s variant_calling.rules --forceall

.. seealso:: :ref:`pipelines` section for more information.

Using **sequana** standalone
------------------------------

An easier way to initialise a pipeline, is to use **sequana** executable. For
instance for the variant calling::

    sequana --pipeline variant_calling

This will automatically download the pipeline, config file and update the latter
as much as possible.

.. seealso:: :ref:`applications` section


Using **Sequanix** standalone
---------------------------------

An even easier way is to use our graphical interface named **Sequanix**. A
snapshot can be found in the :ref:`sequanix` section and a tutorial in
:ref:`tutorial_sequanix`.



**Sequana** Reports
=====================


Pipelines and standalone make use of internal reporting. Since there are part of
the **Sequana** library, they can also be used with your own code. For instance,
if you have a BAM file, you can use the following code to create a basic
report::

    from sequana import BAM, sequana_data
    from sequana.modules_report.bamqc import BAMQCModule
    filename sequana_data("test.bam", "testing")

    r = BAMQCModule(filename, "bam.html")

that results can be shown in `bam.html <_static/bam.html>`_

