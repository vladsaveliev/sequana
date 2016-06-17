# Import -----------------------------------------------------------------------

import os
from sequana.report_mapping import MappingReport
from sequana import bedtools, sequana_data

# Test -------------------------------------------------------------------------

def test_report():
    mydata = bedtools.Genomecov(sequana_data("test_bedcov.bed"))
    mydata.running_median(n=3, circular=False)
    mydata.coverage_scaling()
    mydata.compute_zscore()
    r = MappingReport()
    r.set_data(mydata)
    r.create_report()
