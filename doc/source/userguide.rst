User Guide
############


Snakemake tools
================

We provide so utilities in the :mod:`sequana.snakemake` to be used within
Snakefile. Although Python code can be added the print command output may
interfer with the Snakefile interpretation. The :func:`sequana.snakemake.message` is provided to prevent that issue (see example below).


Another feature of Sequana is to provide a set of modules that are independent
and generic enough. For instance, Snakemake has a nice option to create the dag
file out of the snakefile. However, one has to then use **dot** command and
possibly other commands. In Sequana, we also have a dot parser to annotate the
dot (useful in HTML document for clickable image). The module that takes care of
that is called ... dag and can be included in Your Snakefile as follows::

    from sequana import snakemake as sm

    sm.message("Include dag rule in this pipeline" )

    include: sm.rules['dag']

Save this example in a file called Snakefile and type::

    snakemake


