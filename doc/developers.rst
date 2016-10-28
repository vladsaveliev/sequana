.. _developers:

Developer guide
################

In this section we first look at how to include a new module (snakemake rule) in
Sequana. Then, we will create a new pipeline that uses that single rule.

The rule simply counts the number of reads in a fastq file.
The pipeline will contains that unique rule only. 


How to write a new module
==============================

The :term:`Module` in Sequana parlance is a directory with two files: one that
contains a Snakemake file (also known as :term:`Snakefile`) and a README file.
The Snakefile may be a simple snakemake rule or a set of them (a pipeline). 

Find a valid name
-------------------

All rules and pipelines must have a unique name in Sequana.
We can quickly check that a name is not used in Sequana using:

.. doctest::

    >>> import sequana
    >>> "count" in sequana.modules.keys()
    False

So, let us name it **count**


Create a Snakefile rule
-------------------------
A possible code that implements the **count** rule is the following Snakefile:

.. literalinclude:: Snakefile
    :lines: 3-20
    :linenos:
    :language: python
    :emphasize-lines: 4


This is not a tutorial on Snakemake but let us quickly explain this Snakefile.
The first two lines use **Sequana** library to provide the *filename* as a test file.


Then, the rule itself is defined on line 3 where we define the rule name: **count**. We 
then provide on line 5 and 6 the expected input and output filenames. On line 7 onwards, we
define the actual function that counts the number of reads and save the results
in a json file.

You can now execute the Snakefile just to check that this rule works as expected::

    snakemake -s Snakefile -f

You can check that the file *count.txt* exists.

.. note:: The option `-f ` forces snakemake to run the rules (even though it was 
    already computed earlier).


Store the rule in a Sequana module
-------------------------------------

We now store this Snakefile in the proper place. 
All modules are placed either in *./sequana/pipelines* or in
*./sequana/rules* directory. The tree structure looks like::

    .
    ├── config.yaml
    ├── rules
    │   ├── count
    │   │   ├── count.rules
    │   │   ├── README.rst
    ├── pipelines
    │   ├── count_pipeline
    │   │   ├── count_pipeline.rules
    │   │   ├── README.rst


We have created a *count* directory in *./rules* and put the Snakefile in it
(named *count.rules*).

A few comments: (1) there is a unique file named *config.yaml* at the root, (2)
directory names must match the rule filename contain in it (3) 
the Snakefiles all end in *.rules*, and (4) a *README.rst* must be present in all
directories. 

The README file in the rules can be empty. However, the README in the
pipelines's directory is used in the documentation and automatically
parsed. See :ref:`readme` section for further details.

The count rules is now part of the library, which can be checked using the same
code as before:

::

    >>> import sequana
    >>> "count" in sequana.modules.keys()
    True


Convention to design a rule 
======================================


Use variables
-------------------------------------

Consider this Snakefile::

    rule bedtools_genomecov:
        input:
            __bedtools_genomecov__input
        output:
            __bedtools_genomecov__output
        params:
            options = config["bedtools"]["options"]
        shell:
            bedtools genomecov {params.options} -ibam {input} > {output}

We tend to not hard-code any filename. So the input and output are actually
variables. The variable names being the name of the rule with leading and
trailing doubled underscores followed by the string *input* or *output*.

The advantage is that those variables can be defined in the Snakefile before the
rule but also overwritten within a pipeline as discussed in :ref:`dev_pipeline`

Use the config file
---------------------------------

We encourage developers to NOT set any parameter in the params section of the Snakefile.
Instead , put all parameters required inside the *config.yaml* file.
Since each rule has a unique name, we simply add a section with the rule name. For 
instance::

    bedtools_genomecov:
        options: ''

This is a YAML formatted file. Note that there is no information here. However, 
one may provide any parameters understood by the rule (here *bedtools genomecov*
application) in the *options* field. 

We encourage developers to put as few parameters as possible inside the config.
First to not confuse users and second because software changes with time. Hard
coded parameter may break the pipeline. However having the *options* field
allows users to use any parameters.


.. seealso:: Sequana contains many pipelines that can be used as examples.
    See `github repo <https://github.com/sequana/sequana/tree/master/sequana/pipelines>`_


Add documentation
------------------------

In **sequana**, we provide a sphinx extension to include the inline
documentation of a rule::

    .. snakemakerule:: rule_name

This searches for the rule docstring, and includes it in your documentation. The
docstring should be uniformised across all rules and pipelines. Here is our
current convention::

    rule cutadapt:
        """Cutadapt (adapter removal)

        Some details about the tool on what is does is more than welcome.

        Required input:
            - __cutadapt__input_fastq

        Required output:
            - __cutadapt__output

        Required parameters:
            - __cutadapt__fwd: forward adapters as a file, or string
            - __cutadapt__rev: reverse adapters as a file, or string

        Required configuration:
            .. code-block:: yaml

                cutadapt:
                    fwd: "%(adapter_fwd)s"
                    rev: "%(adapter_rev)s"

        References:
            a url link here or a link to a publication.
    """
    input:
        fastq = __cutadapt__input_fastq
    output:
        fastq = __cutadapt__output
    params:
        fwd = config['cutadapt']["fwd"],
        rev = config['cutadapt']["rev"]
    run:
        cmd = "cutadapt -o {output.fastq[0]} -p {output.fastq[1]} "
        cmd += " -g %s -G %s" % (params.fwd, params.rev)
        cmd += "{input.fastq[0]} {input.fastq[1]}"
        shell(cmd)


.. _dev_pipeline:

How to write new pipelines
================================

Use sequana.snaketools
-------------------------

Assuming there will be a config file named *config.yaml*, the pipeline should be
written as follows:

.. code-block:: snakemake

    import sequana
    from sequana import snaketools as sm
    sm.init("counter", globals())        # see later for explanation

    configfile: "config.yaml"

    # include all relevant rules
    include: sm.modules['count']        # if included in sequana/rules


    # must be defined as the final rule
    rule pipeline_count:
        input:
            "count.txt"


.. _readme:

The pipeline README file
----------------------------

This is a WIP in progress but here is an example for the previous pipeline::

    :Overview: Counts the reads in a fastq file
    :Input: FastQ raw data file
    :Output: 
        - count.txt

    Usage
    ~~~~~~~

    ::

        sequana init pipeline_count 
        snakemake -s pipeline_count.rules -f 

.. note:: the README uses Restructured syntax (not markdown)




requirements file
-----------------------

You can add a **requirements.txt** file with a list of executable or libraries
required by the pipeline. Executable are searched in your PATH and libraries in
PYTHONPATH (conda).

Where running a pipeline, if the requirements are not fulfilled, an error is
raised before starting the analysis.




