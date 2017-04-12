from sequana.pacbio import BAMPacbio
from sequana import sequana_data
from easydev import TempFile

def test_pacbio():
    b = BAMPacbio(sequana_data("test_pacbio_subreads.bam"))
    assert len(b) == 8
    b.df
    assert b.nb_pass[1] == 8

    with TempFile() as fh:
        b.stride(fh.name, stride=2)

    with TempFile() as fh:
        b.filter_length(fh.name, threshold_min=500)

    b.hist_snr()
    b.hist_ZMW_subreads()

    b.hist_len()
    b.hist_GC()
    b.plot_GC_read_len()
    
