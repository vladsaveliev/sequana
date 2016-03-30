from sequana.bamtools import BAM, Alignment, SAMFlags, BAMReport
from sequana import sequana_data
from easydev import TempFile



def test_bam():
    s = BAM(sequana_data("test.bam"))
    assert len(s) == 1000

    assert len(list(s.iter_unmapped_reads())) == 2
    s.reset()
    assert len(list(s.iter_mapped_reads())) == 998
    s.reset()

    assert s.get_read_names()

    s.get_stats()

    with TempFile(suffix='.png') as fh:
        s.plot_bar_flags(filename=fh.name, logy=True)
        s.plot_bar_flags(filename=fh.name)
    

    with TempFile(suffix='.png') as fh:
        s.plot_bar_mapq(filename=fh.name)
  

def test_alignment(): 
    s = BAM(sequana_data("test.bam"))
    a = Alignment(next(s))
    a.as_dict()


def test_samflags():

    sf = SAMFlags(4095)
    print(sf)
    sf.get_meaning()
    sf.get_flags()


def test_bamreport():
    b = BAM(sequana_data("test.bam"))
    r = BAMReport()
    r.set_data(b)
    r.create_report()
