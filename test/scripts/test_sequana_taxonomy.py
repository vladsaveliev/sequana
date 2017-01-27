from sequana.scripts import taxonomy
from sequana import sequana_data
import pytest

prog = "sequana_taxonomy"

import os

if "TRAVIS_PYTHON_VERSION" not in os.environ:
    @pytest.fixture
    def krakendb():
        # todo
        try:
            taxonomy.main(["taxonomy", '--download', 'toydb'])
        except SystemExit:
            pass


    def test_analysis(krakendb):
        file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
        file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz")

        # Test that database must be provided
        try:
            df = taxonomy.main([prog, '--file1', file1])
            assert False
        except:
            assert True

        from tempfile import TemporaryDirectory
        directory = TemporaryDirectory()
        df = taxonomy.main([prog, '--file1', file1, "--database", "toydb",
            "--file2", file2, "--verbose", "--output-directory", directory.name])
        from sequana import logger
        logger.info(directory.name)

    def test_help():
        try:
            taxonomy.main([prog, '--help', '1>/tmp/out', '2>/tmp/err'])
            assert False
        except SystemExit:
            pass
        else:
            raise Exception

def _test_wrong_db():
    try:
        df = taxonomy.main([prog, "--database", "dummy"])
        assert False
    except:
        assert True
