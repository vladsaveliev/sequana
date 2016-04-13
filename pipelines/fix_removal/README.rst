Usage
=======


::

    sequana --init  phix_removal

Edit the config.yaml to update the pattern for the input files (first line)


Run as follows:: 

    snakemake -j 2 -p --stats stats.txt

Or type::

    sh sequana.sh
