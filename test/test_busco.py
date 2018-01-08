from sequana.busco import BuscoConfig, BuscoDownload
from sequana import sequana_data
from easydev import TempFile

def test_busco_config():
    bc = BuscoConfig("species", outpath="test", sample_name="test",
            conda_bin_path="test", tmp_path="test",
            Rscript_bin_path=None)
    with TempFile() as fh:
        bc.save_config_file(fh.name)

def test_busco_download():
    bd = BuscoDownload()
    bd.filenames = ['proteobacteria_odb9']
    bd.download()
