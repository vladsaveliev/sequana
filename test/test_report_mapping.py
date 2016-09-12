import os

from sequana.reporting.report_mapping import MappingReport
from sequana.reporting.report_chromosome import ChromosomeMappingReport
from sequana import bedtools, sequana_data


def test_report():
    mydata = bedtools.GenomeCov(sequana_data("test_bedcov.bed"))
    r = MappingReport()
    r.set_data(mydata)
    r.create_report()
    for chrom in mydata:
        chrom.running_median(n=501, circular=False)
        chrom.compute_zscore()
        r = ChromosomeMappingReport(chrom)
        r.create_report()
