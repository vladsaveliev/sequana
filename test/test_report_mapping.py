# Import -----------------------------------------------------------------------

import os
from sequana.report_mapping import MappingReport
from sequana import bedtools
from . import data
pathdata = data.__path__[0]

# Test -------------------------------------------------------------------------

def test_report():
    mydata = bedtools.Genomecov(pathdata + os.sep + "test.bed")
    mydata.running_median(n=3, circular=False)
    mydata.coverage_scaling()
    mydata.compute_zscore()
    r = MappingReport()
    r.set_data(mydata)
    r.create_report()
