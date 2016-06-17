"""
BAM module example
====================

Histogram MAPQ
"""
from pylab import *

#################################################
# ignore warning message for now
from sequana import BAM, sequana_data


#####################################################
#
from sequana import BAM, sequana_data
b = BAM(sequana_data('test.bam', "testing"))
b.plot_bar_mapq()
#show()
