import os

from sequana import bedtools, sequana_data
from sequana.tools import genbank_features_parser
from easydev import TempFile


def test_threshold():
    t = bedtools.DoubleThresholds(-5,5)
    assert t.low == -5
    assert t.high == 5
    assert t.low2 == -2.5
    t = bedtools.DoubleThresholds(-4, 3)
    assert t.low == -4
    assert t.high == 3

    t = bedtools.DoubleThresholds(-8,8)
    t.ldtr = 0.25
    t.hdtr = 0.25
    assert t.low2 == -2
    assert t.high2 == 2
    print(t)


    t.ldtr = 0.5
    t.hdtr = 0.5
    t.low = -3
    t.high = 3
    assert t.low2 == -1.5
    assert t.high2 == 1.5

    try:
        t = bedtools.DoubleThresholds(3, 4)
        assert False
    except:
        assert True
    try:
        t = bedtools.DoubleThresholds(3, -4)
        assert False
    except:
        assert True

def test_genomecov():

    filename = sequana_data('JB409847.bed')
    bed = bedtools.GenomeCov(filename, sequana_data('JB409847.gbk'))

    # a getter for the first chromosome
    bed[0]

    # This requires to call other method before
    for chrom in bed:
        chrom.moving_average(n=501)
        chrom.running_median(n=501, circular=True)
        chrom.running_median(n=501, circular=False)

        chrom.compute_zscore()
        roi = chrom.get_roi()
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
    with TempFile(suffix='.csv') as fh:
        bed.to_csv(fh.name)
        bed2 = bedtools.GenomeCov(fh.name, sequana_data('JB409847.gbk'))

def test_gc_content():
    bed = sequana_data('JB409847.bed')
    fasta = sequana_data('JB409847.fasta')
    cov = bedtools.GenomeCov(bed)
    cov.compute_gc_content(fasta)
    cov.get_stats()
    ch = cov[0]
    ch.moving_average(4001, circular=True)
    ch.running_median(4001,circular=True)
    ch.compute_zscore()

    ch.get_evenness()
    ch.get_cv()
    assert ch.get_centralness() > 0.84 and ch.get_centralness()<0.85
    ch.plot_gc_vs_coverage()

    from easydev import TempFile
    with TempFile() as fh:
        ch.to_csv(fh.name)


    ch.get_max_gc_correlation(fasta)
