from sequana import KronaMerger, KrakenResults, KrakenTaxon, KrakenBuilder
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


@attr('skip')
def test_kraken_taxon():

    from easydev import execute
    import tempfile

    tmpdir = tempfile.mkdtemp()
    execute("git clone http://github.com/sequana/data %s" % tmpdir)


    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "testing")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "testing")

    kt = KrakenTaxon([file1, file2], database="%s/kraken_toydb" % tmpdir)
    kt.run()

    kt = KrakenTaxon(file2, database="%s/kraken_toydb" % tmpdir)
    kt.run()

    #Clean the kraken db
    import shutil
    shutil.rmtree(tmpdir)


def _test_kraken_builder():


    k = KrakenBuilder(download_taxon=True)
    k.run(['viroids'])

