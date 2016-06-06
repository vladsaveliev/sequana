How to write pipelines
==============================

New pipelines must be added in ./sequana/pipelines directory. 
In sequana, all rules and pipelines must be uniquely named. So first identify 
a name that is not already used. For this example, let us use **newpipe**.

By convention, all pipelines and rules ends with **.rules** extension and a
unique config file named **config.yaml** is stored in
sequana/pipelines directory.

The new pipeline must be stored in with this tree structure::

    .
    ├── config.yaml
    ├── pipelines
    │   ├── newpipe
    │   │   ├── newpipe.rules
    │   │   ├── README.rst



The content of **newpipe.rules** must specify the name of the pipeline itself,
which is going to be used in the **dag.rules** ::

    __snakefile__ = "newpipe.rules"


.. seealso:: Sequana contains many pipelines that can be used as examples.


How to write rules
========================


.. seealso:: Sequana contains many rules that can be used as examples.

