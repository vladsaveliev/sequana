"""
Coverage module example
=============================


"""
#################################################
#
from sequana import GenomeCov
from sequana import sequana_data
bedfile = sequana_data("virus.bed", "data")
reference = sequana_data("pKM2-eGFP.fa", "data")
################################################
# Reading input BED file
gc = GenomeCov(bedfile)

##############################################
# Reading reference and computing GC content
w_gc = 100
gc.compute_gc_content(reference, w_gc)


##########################################################
# Select a chromosome (first and only one in that case
chrom = gc[0]

##########################################################
# Compute running median and zscore
chrom.running_median(n=4001, circular=True)
chrom.compute_zscore()

######################################"##################
# Plotting
chrom.plot_coverage()


