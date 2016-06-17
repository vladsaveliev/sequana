from sequana import KronaMerger 
from sequana import sequana_data
from easydev import TempFile, execute
from nose.plugins.attrib import attr

def test_krona_merger():

    k1 = KronaMerger(sequana_data("test_krona_k1.tsv", "testing"))
    k2 = KronaMerger(sequana_data("test_krona_k2.tsv", "testing"))
    k1 += k2
 
    with TempFile(suffix='.tsv') as fh:
        df = k1.to_tsv(fh.name) 
    assert all(df['count'] == [14043,591,184,132])
    assert k1['Bacteria\tProteobacteria\tspecies1\n'] == 14043


