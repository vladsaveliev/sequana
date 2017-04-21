.. _developers:

Developer guide
################

In this section we first look at how to include a new module (snakemake rule) in
**Sequana**. Then, we will create a new pipeline that uses that single rule.

The rule simply counts the number of reads in a fastq file.
The pipeline will only contains that unique rule. 

In the remaining sections, we will explain our choice concerning the continuous
integration (section :ref:`pytest` ) and how to add check that new code do not
introduce bugs. In the :ref:`module_reports` section, we explain how to create
new component in the HTML module reports.


How to write a new module
==============================

A :term:`Module` (in Sequana parlance) is a directory with a set of files:
a Snakemake file (also known as :term:`Snakefile`), a README file for the
documentation and a configuration file (optional).
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


Then, the rule itself is defined on line 4 where we define the rule named: **count**. We 
then provide on line 5 and 6 the expected input and output filenames. On line 7 onwards, we
define the actual function that counts the number of reads and save the results
in a TXT file.

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

.. note:: The is a big advantage of designing rules with variables only:
    rules can be re-used in any pipelines without changing the rule itself; only 
    pipelines will be different.

Use a config file
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


Add documentation in the rule
----------------------------------

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

There are many rules already available in **Sequana**. You can easily add rules
as follows::

    from sequana import snaketools as sm
    include: sm.module['rulegraph']


Use sequana.snaketools
-------------------------

Assuming there will be a config file named *config.yaml*, the pipeline should be
written as follows:

.. code-block:: python

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


Documenting the configuration file
---------------------------------------

The configuration should be in YAML format. You should comment top-level
sections corresponding to a rule as follows::


    #######################################################
    #
    # A block comment in docstring format
    #
    # This means a # character followed by a space and then 
    # the docstring. The first line made of ##### will be removed
    # and is used to make the documentation clear. No spaces
    # before the section (count:) here below.
    #
    count:  # you can add an overview
        -item1: 1 # you can add comment for an item
        -item1: 2 # you can add comment for an item

If valid, the block comment is interpreted and a tooltip will appear in
**Sequanix**.


Further coding conventions
-------------------------------

To print debugging information, warnings or more generally information, please
do not use the print() function but the logger::

    from sequana import logger
    logger.debug("test")
    logger.info("test")
    logger.warning("test")
    logger.error("test")
    logger.critical("test")


.. _pytest:

Testing with pytest
===============================

We moved from nosetests to pytest. This framwork is slightly more flexible but
the main reason to move was to be able to test Qt application. 

In order to run the test locally, you will need to use::

    pip install pytest pytest-cov pytest-qt pytest-xdist pytest-mock pytest-timeout


Then, you can type for instance::

    pytest -v --durations=10  test/ --cov=sequana --cov-report term --timeout 300 -n 2

Here, -n 2 requires two CPUs to run the test. The option durations=10 means
*show the 10 longest tests*. 


If you want to test a single file (e.g. test_pacbio)::

    cd test
    pytest test_pacbio.py --cov sequana.pacbio --cov-report term-missing



.. _module_reports:

Module reports
======================

Sequana pipelines generate HTML reports. Those reports are created with
the **module reports** stored in ./sequana/modules_report directory.

A module report creates one HTML page starting from a dataset generated
by Sequana, or a known data structure. All modules reports inherit from
:class:`~sequana.modules_report.base_module.SequanaBaseModule` as shown
hereafter. This class provides convenient methods to create the final HTML,
which takes care of copying CSS and Javascript libraries.


To explain how to write a new module report, let us consider a simple example. 
We design here below a working example of a module report that takes as input a
Pandas dataframe (a Pandas series made of a random normal distribution to be precise).
The module report then creates an HTML page with two sections: a dynamic sortable table
and a section with an embedded image. Each section is made of a dictionary that
contains 3 keys:

- name: the HTML section name
- anchor: the ID HTML of the section
- content: a valid HTML code

First, you need to import the base class. Here we also import a convenient
object called DataTable that will be used to created sortable table in HTML
using javascript behind the scene. 

.. code-block:: python

    from sequana.modules_report.base_module import SequanaBaseModule
    from sequana.utils.datatables_js import DataTable

Then, we define a new class called **MyModule** as follows::

    class MyModule(SequanaBaseModule):

followed by a constructor


.. literalinclude:: module_example.py
    :language: python
    :pyobject: MyModule.__init__

This constructor stores the input argument (**df**) and computes some new data
stored in the :attr:`summary` attribute. Here this computation is fast but in
a real case example where computation may takes time, the
computation  should be performed outside of the module. We then store a title
in the :attr:`title` attribute. Finally two methods are called. The first one
creates the HTML sections (*create_report_method*); the second one
(*create_html*) is inherited from SequanaBaseModule.

The first method is defined as follows:

.. literalinclude:: module_example.py
    :language: python
    :linenos:
    :pyobject: MyModule.create_report_content
    :emphasize-lines: 2

Here, the method :meth:`create_report_content` may be named as you wish but must
define and fill the :attr:`sections` list (empty list is possible) with a set
of HTML sections. In this example, we call two methods (add_table and add_image)
that adds two HTML sections in the list. You may have as many **add_** methods.

First, let us look at the :meth:`add_table`. It creates an HTML section made of
a dynamic HTML table based on the DataTable class. This class takes as input a
Pandas DataFrame.

.. literalinclude:: module_example.py
    :language: python
    :linenos:
    :pyobject: MyModule.add_table


Here, we first get some data (line 2) in the form of a Pandas time Series. We
rename the column on line 3. This is a dataframe and the DataTable class takes
as input Pandas dataframe that are then converted into flexible HTML table.

One nice feature about the DataTable is that we can add HTML links (URL) in
a specific column of the data frame (line 3) and then link an existing column with this
new URL column. This happens on line 11. The final HTML table
will not show the URL column but the data column will be made of clickable
cells.

The creation of the data table itself happens on line 5 to line 11 and line
12-14. There are two steps here: the creation of the HTML table itself (line 13)
and the Javascript itself (line 12). 

Once we have the HTML data, we can add it into the sections on line 16-19.

The second section is an HTML section with an image. It may be 
included with a standard approach (using the **img** tag) but one can also use the
:meth:`~sequana.modules_report.SequanaBaseModule.create_embedded_png` method.


.. literalinclude:: module_example.py
    :language: python
    :linenos:
    :pyobject: MyModule.add_image

Here is the full working example: 

.. literalinclude:: module_example.py
    :language: python
    :linenos:


When using this module, one creates an HTML page called **mytest.html**. An
instance of the page is available here:  `report_example.html <_static/report_example.html>`_

