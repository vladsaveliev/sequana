import os

from sequana import bedtools, sequana_data
from easydev import TempFile


def test_genomecov():

    filename = sequana_data("test_bedcov.bed", "testing")

    mydata = bedtools.GenomeCov(filename)
  
    # This requires to call other method before
    for chrom in mydata:
        chrom.moving_average(n=501)
        chrom.running_median(n=501, circular=True)
        chrom.running_median(n=501, circular=False)

        chrom.compute_zscore()
        chrom.get_low_coverage()
        high_cov = chrom.get_high_coverage()
        high_cov.merge_region(3)
        with TempFile(suffix='.png') as fh:
            chrom.plot_coverage(filename=fh.name)
        with TempFile(suffix='.png') as fh:
            chrom.plot_hist_zscore(filename=fh.name)
        with TempFile(suffix='.png') as fh:
            chrom.plot_hist_normalized_coverage(filename=fh.name)

        len(chrom)
        print(chrom)
