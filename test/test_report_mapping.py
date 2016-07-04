import os

from sequana.report_mapping import MappingReport
from sequana.report_chromosome import ChromosomeMappingReport
from sequana import bedtools, sequana_data


def test_report():
    mydata = bedtools.GenomeCov(sequana_data("test_bedcov.bed"))
    r = MappingReport()
    r.set_data(mydata)
    r.create_report()
    for chrom in mydata:
        chrom.running_median(n=3, circular=False)
        chrom.coverage_scaling()
        chrom.compute_zscore()
        r = ChromosomeMappingReport(chrom=chrom.chrom_name)
        r.set_data(chrom)
        r.create_report()
