from sequana.scripts import compressor
from sequana import sequana_data, FastQ
import os
import tempfile
from easydev import TempFile
import easydev
import shutil
import argparse
from unittest.mock import MagicMock, patch

prog = "sequana_compressor"

def test_compressor_args(mocker):

    argparse_mock = MagicMock()
    with patch('argparse.ArgumentParser._print_message', argparse_mock):
        try:
            compressor.main(["sequana_compressor", '--version'])
            assert False
        except SystemExit:
            assert True
        try:
            compressor.main(["sequana_compressor"])
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

    try:
        compressor.main([prog, "--jobs", "25"])
        assert False
    except ValueError:
        assert True


    class Node():
        node = "tars-submit"
    return Node()
    with patch('platform.uname', Node()):
        try:
            compressor.main([prog, "--source", "fastq.gz", "--target", "fastq.bz2"])
            assert False
        except ValueError:
            assert True


def test_compressor_bad_extension():
    try:
        compressor.main([prog, "--source", "fastq.txt", "--target", "fastq.txt"])
        assert False
    except ValueError:
        pass
    else:
        raise Exception


def test_compressor_running():
    # Here we test gz -> bz2 -> gz -> dsrc -> gz using recursive or not
    # get a fastq.gz in a temp file and process it
    tempdir = tempfile.TemporaryDirectory()
    filename = sequana_data("test.fastq.gz")
    shutil.copy(filename, tempdir.name )
    cwd = os.path.abspath(os.curdir)
    os.chdir(tempdir.name)

    # We concert gz -> bz2 -> gz and must get the exact same
    # However, since the compression is not deterministic, we should compare the
    # content of the uncompressed file (input and output)
    try:
        # seems to fail on travis with a subprocess issue
        # https://travis-ci.org/sequana/sequana/builds/162466158
        compressor.main([prog, "--source", "fastq.gz", "--target", "fastq.bz2", "--quiet" ])
        compressor.main([prog, "--source", "fastq.bz2", "--target", "fastq.gz", "--recursive", "--quiet"])
        compressor.main([prog, "--source", "fastq.gz", "--target", "fastq.dsrc", "--recursive", "--quiet"])
        compressor.main([prog, "--source", "fastq.dsrc", "--target", "fastq.gz", "--quiet"])
    except Exception as err:
        raise Exception(err)
    finally:
        os.chdir(cwd)

    f1 = FastQ(filename)
    f2 = FastQ(tempdir.name + os.sep + os.path.basename(filename))
    assert f1 == f2


def test_compressor_dsrc():
    """Test dsrc codecs to gz and bz2"""
    # Create a temporary directory and chdir in it:
    tempdir = tempfile.TemporaryDirectory()
    filename = sequana_data("test.fastq.gz")
    shutil.copy(filename, tempdir.name )
    cwd = os.path.abspath(os.curdir)
    os.chdir(tempdir.name)

    # We concert gz -> dsrc -> bz2 -> dsrc -> gz and must get the exact same
    # However, since the compression is not deterministic, we should compare the
    # content of the uncompressed file (input and output)
    try:
        compressor.main([prog, "--source", "fastq.gz", "--target", "fastq.dsrc", "--quiet" ])
        compressor.main([prog, "--source", "fastq.dsrc", "--target", "fastq.bz2", "--quiet"])
        compressor.main([prog, "--source", "fastq.bz2", "--target", "fastq.dsrc", "--quiet"])
        compressor.main([prog, "--source", "fastq.dsrc", "--target", "fastq.gz", "--quiet"])
    except Exception as err:
        raise Exception(err)
    finally:
        os.chdir(cwd)

    f1 = FastQ(filename)
    f2 = FastQ(tempdir.name + os.sep + os.path.basename(filename))
    assert f1 == f2












