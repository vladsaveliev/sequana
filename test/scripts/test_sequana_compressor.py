from sequana.scripts import compressor
from sequana import sequana_data, FastQ
import os
import tempfile
from easydev import TempFile
import easydev
import shutil
import argparse
from unittest.mock import MagicMock, patch


def test_compressor_args():
    prog = "sequana_compressor"

    argparse_mock = MagicMock()
    with patch('argparse.ArgumentParser._print_message', argparse_mock):
        try:
            compressor.main(["sequana_compressor", '--version'])
            assert False
        except SystemExit:
            assert True
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


def test_compressor_running():
    prog = "sequana_compressor"

    # get a fastq.gz in a temp file and process it
    tempdir = tempfile.TemporaryDirectory()
    filename = sequana_data("test.fastq.gz")
    checksum1 = easydev.md5(filename)
    shutil.copy(filename, tempdir.name )


    cwd = os.path.abspath(os.curdir)
    os.chdir(tempdir.name)

    try:
        # seems to fail on travis with a subprocess issue
        # https://travis-ci.org/sequana/sequana/builds/162466158
        compressor.main([prog, "--source", "fastq.gz", "--target", "fastq.bz2", "--quiet" ])
        compressor.main([prog, "--source", "fastq.bz2", "--target", "fastq.gz", "--recursive", "--quiet"])
        compressor.main([prog, "--source", "fastq.gz", "--target", "fastq.dsrc", "--recursive", "--quiet"])
        compressor.main([prog, "--source", "fastq.dsrc", "--target", "fastq.gz", "--quiet"])
    except err:
        raise Exception(err)
    finally:
        os.chdir(cwd)

    f1 = FastQ(filename)
    f2 = FastQ(tempdir.name + os.sep + os.path.basename(filename))
    assert f1 == f2





