# Import -----------------------------------------------------------------------

import os
from sequana import bedtools, sequana_data
from easydev import TempFile

# Test -------------------------------------------------------------------------

def test_genomecov():

    filename = sequana_data("test_bedcov.bed", "testing")

    mydata = bedtools.GenomeCov(filename)
  
    # This requires to call other method before
    for chrom in mydata:
        chrom.coverage_scaling()
        chrom.compute_zscore()

        chrom.moving_average(n=101)
        chrom.running_median(n=101, circular=True)
        chrom.running_median(n=101, circular=False)
        chrom.coverage_scaling()

        chrom.compute_zscore()
        chrom.get_low_coverage()
        high_cov = chrom.get_high_coverage()
        high_cov.merge_region(101)
        with TempFile(suffix='.png') as fh:
            chrom.plot_coverage(filename=fh.name)
        with TempFile(suffix='.png') as fh:
            chrom.plot_hist(filename=fh.name)

        len(chrom)
        print(chrom)
