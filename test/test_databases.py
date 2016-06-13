from sequana import databases
import os

def test_database_download():

    d = databases.ENADownload()
    d.download_viroid()
    import glob
    for this in glob.glob('Viroid/*'):
        os.remove(this)
    
    os.rmdir("Viroid")


    d.download_list()
    # cleanup
    for key in d._metadata.keys():
        res = d._metadata[key]
        filename = res[0]
        os.remove(filename)


    d.download_accession('AF439431')
    os.remove("Custom/ENA_AF439431.fa")
    os.rmdir("Custom")

