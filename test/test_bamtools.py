from sequana.bamtools import BAM, Alignment, SAMFlags
from sequana.modules_report.bamqc import BAMQCModule
from sequana import sequana_data
from easydev import TempFile


datatest = sequana_data("test.bam", "testing")


def test_bam(tmpdir):

    s = BAM(datatest)
    assert len(s) == 1000
    assert s.is_sorted is False

    assert len(list(s.iter_unmapped_reads())) == 2
    s.reset()
    assert len(list(s.iter_mapped_reads())) == 998
    s.reset()

    # call this here before other computations on purpose
    with TempFile(suffix=".json") as fh:
        s.bam_analysis_to_json(fh.name)

    assert s.get_read_names()
    s.get_mapped_read_length()

    s.get_stats()
    s.get_full_stats_as_df()

    with TempFile(suffix='.png') as fh:
        s.plot_bar_flags(filename=fh.name, logy=True)
        s.plot_bar_flags(filename=fh.name)

    with TempFile(suffix='.png') as fh:
        s.plot_bar_mapq(filename=fh.name)

    with TempFile() as fh:
        s.to_fastq(fh.name)
        from sequana import FastQ
        ff = FastQ(fh.name)
        len(ff) == len(s)

    s.get_gc_content()
    s.get_length_count()
    s.plot_gc_content()
    try:
        s.plot_gc_content(bins=[1,2,10])
        assert False
    except:
        assert True


def test_bam_others():
    b = BAM(sequana_data("measles.fa.sorted.bam"))
    assert len(b) == 2998

    # plot_reaqd_length and data
    X, Y = b._get_read_length()
    assert sum(Y) == 2623
    b.plot_read_length()
    b.plot_acgt_content()
    b.hist_coverage()
    b.plot_coverage()
    b.boxplot_qualities()
    b.plot_indel_dist()    

def test_alignment():
    s = BAM(datatest)
    # no need to call reset but does not harm and reminds us that it shoudl be
    # used in general to make sure we start at the beginning of the iterator.
    s.reset()
    a = Alignment(next(s))
    a.as_dict()


def test_samflags():
    sf = SAMFlags(4095)
    print(sf)
    sf.get_meaning()
    sf.get_flags()


def test_bamreport(tmpdir):
    directory = tmpdir.mkdir("bam")
    from sequana.utils import config
    config.output_dir = directory.__str__()
    r = BAMQCModule(datatest, "bam.html")
