import os

from sequana import bedtools, sequana_data
from sequana.tools import genbank_features_parser
from easydev import TempFile


def test_genomecov():

    filename = sequana_data("test_bedcov.bed", "testing")
    features = genbank_features_parser(sequana_data("test_snpeff_ref.gb"))

    mydata = bedtools.GenomeCov(filename)
  
    # This requires to call other method before
    for chrom in mydata:
        chrom.moving_average(n=501)
        chrom.running_median(n=501, circular=True)
        chrom.running_median(n=501, circular=False)

        chrom.compute_zscore()
        roi = chrom.get_roi(features=features)
        with TempFile(suffix='.png') as fh:
            chrom.plot_coverage(filename=fh.name)
        with TempFile(suffix='.png') as fh:
            chrom.plot_hist_zscore(filename=fh.name)
        with TempFile(suffix='.png') as fh:
            chrom.plot_hist_normalized_coverage(filename=fh.name)

        len(chrom)
        print(chrom)
        chrom.get_size()
        chrom.get_mean_cov()
        chrom.get_var_coef()
