from sequana.samtools import SAM

import os
from . import data
pathdata = data.__path__[0]



def test_sam():
    s = SAM(pathdata + os.sep + "test.sam")
    assert len(s.get_read_names()) == 432
    s.plot_mapq_distribution()
