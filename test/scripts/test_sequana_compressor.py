from sequana.scripts import compressor
from nose.plugins.attrib import attr
from sequana import sequana_data
import os

#@attr("skip")
def test_compressor():
    prog = "sequana_compressor"
    try:
        compressor.main(["sequana_compressor", '--version'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception

    try:
        compressor.main([prog, '--help'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception
    try:
        compressor.main([prog])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception


    # get a fastq.gz in a temp file and process it
    from easydev import TempFile
    import shutil
    filename = sequana_data("test.fastq.gz")
    fh = TempFile(suffix=".fastq.gz")
    shutil.copy(filename, fh.name)

    cwd = os.path.abspath(os.curdir)
    os.chdir("/tmp")
    compressor.main([prog, "--source", "fastq.gz", "--target", "fastq.bz2", "--quiet"])
    compressor.main([prog, "--source", "fastq.bz2", "--target", "fastq.gz",
"--recursive", "--quiet"])
    fh.delete()



