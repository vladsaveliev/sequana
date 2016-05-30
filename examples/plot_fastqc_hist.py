"""
Quality histogram a la fastqc
=================================


"""

#####################################
# Description

from sequana import sequana_data
from sequana import FastQC

# Create a FastQC instance
filename  = sequana_data("test.fastq", "testing")
qc = FastQC(filename)

# plot the histogram
qc.boxplot_quality()
