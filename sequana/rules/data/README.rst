:Description: Include fastq_sampling rule (updating the config file on the fly)
    to use a sub-set of the raw data. If not used, the raw data is used as
    expected. Note, however that without this rule, the data is not shown in the
    workflow visual representation (DAG).
:Used by: the pipeline :ref:`pipeline_phix_removal` and children pipelines (e.g.
    :ref:`pipeline_quality`).
:Config sections:
    - fastq_sampling
:Dependencies:
    - sequana.fastq_head executable
