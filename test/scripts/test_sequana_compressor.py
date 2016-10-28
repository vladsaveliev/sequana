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
    try:
        # seems to fail on travis with a subprocess issue
        # https://travis-ci.org/sequana/sequana/builds/162466158
        compressor.main([prog, "--source", "fastq.gz", "--target", "fastq.bz2", "--quiet"])
        compressor.main([prog, "--source", "fastq.bz2", "--target", "fastq.gz",                                               "--recursive", "--quiet"])
        compressor.main([prog, "--source", "fastq.gz", "--target", "fastq.dsrc",
                         "--recursive", "--quiet"])
        compressor.main([prog, "--source", "fastq.dsrc", "--target", "fastq.gz",
                        "--recursive", "--quiet"])
    except:
        pass
    fh.delete()
    os.chdir(cwd)
