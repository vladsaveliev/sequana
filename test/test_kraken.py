from sequana import  KrakenResults, KrakenAnalysis
from sequana import sequana_data
from easydev import TempFile, execute
from nose.plugins.attrib import attr



@attr('skip')
def test_kraken_taxon():

    from easydev import execute
    import tempfile

    tmpdir = tempfile.mkdtemp()
    execute("git clone http://github.com/sequana/data %s" % tmpdir)


    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "testing")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "testing")

    kt = KrakenAnalysis([file1, file2], database="%s/kraken_toydb" % tmpdir)
    kt.run()

    kt = KrakenAnalysis(file2, database="%s/kraken_toydb" % tmpdir)
    kt.run()

    #Clean the kraken db
    import shutil
    shutil.rmtree(tmpdir)


def _test_kraken_builder():

    k = KrakenBuilder(download_taxon=True)
    k.run(['viroids'])

