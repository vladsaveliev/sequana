# Import -----------------------------------------------------------------------

import os
from sequana import bedtools, sequana_data
from easydev import TempFile

# Test -------------------------------------------------------------------------

def test_genomecov():

    filename = sequana_data("test_bedcov.bed", "testing")

    mydata = bedtools.Genomecov(filename)
  
    # This requires to call other method before
    mydata.coverage_scaling()
    mydata.compute_zscore()

    mydata.moving_average(n=101)
    mydata.running_median(n=101, circular=True)
    mydata.running_median(n=101, circular=False)
    mydata.coverage_scaling()

    mydata.compute_zscore()
    mydata.get_low_coverage()
    high_cov = mydata.get_high_coverage()
    high_cov.merge_region(101)
    with TempFile(suffix='.png') as fh:
        mydata.plot_coverage(filename=fh.name)
    with TempFile(suffix='.png') as fh:
        mydata.plot_hist(filename=fh.name)



    len(mydata)
    print(mydata)
