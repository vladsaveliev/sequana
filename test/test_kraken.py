from sequana import  KrakenResults, KrakenAnalysis, KrakenDownload
from sequana import sequana_data, sequana_config_path
from nose.plugins.attrib import attr
import os


@attr('skip')
def test_kraken_taxon():

    def download():
        kd = KrakenDownload()
        kd.download('toydb')
    database = sequana_config_path + os.sep + "kraken_toydb"
    if os.path.exists(database) is False:
        download()

    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz", "data")

    kt = KrakenAnalysis([file1, file2], database=database)
    kt.run()

    kt = KrakenAnalysis(file2, database=database)
    kt.run()



def test_kraken_results():
    test_file = sequana_data("test_kraken.out", "testing")
    k = KrakenResults(test_file, verbose=False)
    df = k.plot(kind='pie')
    print(df)

    k = KrakenResults(test_file, verbose=True)
    df = k.plot(kind='barh')

    df = k.get_taxonomy_biokit(11234)
    assert 11234 in df.index

    df = k.get_taxonomy_biokit("11234")
    assert 11234 in df.index


    df = k.get_taxonomy_biokit([11234])
    assert 11234 in df.index

    df = k.get_taxonomy_biokit(["11234"])
    assert 11234 in df.index


def test_download():
    kd = KrakenDownload()
    kd.download('toydb')



