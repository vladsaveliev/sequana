"""
Kraken module example
=======================

kraken module showing distribution of the most frequent taxons
Please, see :mod:`sequana.kraken` for more information and the
quality_taxon pipeline module or kraken rule.
"""
#This plots a simple taxonomic representation of the output
#of the taxonomic pipeline. A more complete and interactive 
#representatino using krona is available when using the 
#quality_taxon pipeline in Sequana.


##############################################
# test 
from sequana import KrakenContaminant
k = KrakenContaminant("kraken.out", verbose=False)
k.plot(kind='pie')


####################################################
# The input file **kraken.out** is the output of the 
# Kraken tool. It is a ste of rows such as those ones::
#
#    C   HISEQ:426:C5T65ACXX:5:2301:5633:7203    11234   203 0:2 11234:1 0:1 11234:1 0:2 11234:1 0:13 11234:1 0:1 11234:1 0:3 11234:1 0:16 11234:1 0:5 11234:1 0:6 11234:1 0:13 A:31 0:33 11234:1 0:29 11234:1 0:7
#    C   HISEQ:426:C5T65ACXX:5:2301:5815:7120    11234   203 0:4 11234:1 0:12 11234:1 0:22 11234:1 0:1 0 11234:1 0:5 11234:1 0:7 11234:1 0:5 A:31 0:3 11234:1 0:22 11234:1 0:18 11234:1 0:24 11234:1
#
#
# The KrakenContaminant class will read the file, download a taxonomic database
# from EBI, map the taxon found in the **kraken.out** file and figure out the
# lineage. In the example above, only the scientific name is found. In the
# snakefile provided in Sequana, the full pipeline produces a full lineage
# representation using krona tool.
#
# .. seealso:: :ref:`pipelines`

