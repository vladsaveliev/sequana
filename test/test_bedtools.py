# Import -----------------------------------------------------------------------

import os
from sequana import bedtools
from easydev import TempFile
from . import data
pathdata = data.__path__[0]

# Test -------------------------------------------------------------------------

def test_genomecov():
    mydata = bedtools.genomecov(pathdata + os.sep + "test.bed")
    assert mydata.running_median(n=3, circular=True)
    assert mydata.running_median(n=3, circular=False)
    assert mydata.coverage_scaling()
    assert mydata.compute_zscore()
    assert mydata.get_low_coverage()
    assert mydata.get_high_coverage()
    assert mydata.merge_region(mydata.df)
    with TempFile(suffix='.png') as fh:
        mydata.plot_coverage(filename=fh)
    with TempDile(suffix='.png') as fh:
        mydata.plot_hist(filename=fh)
