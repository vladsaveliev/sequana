from sequana.bamtools import BAM

import os
from . import data
pathdata = data.__path__[0]


def test_bam():
    s = BAM(pathdata + os.sep + "test.bam")
    assert len(s) == 432

    assert len(list(s.iter_unmapped_reads())) == 432
    s.reset()
    assert len(list(s.iter_mapped_reads())) == 0
    s.reset()

    assert s.get_read_names()[0] == 'HISEQ:426:C5T65ACXX:5:2302:1943:2127'

    s.get_stats()
