from sequana.scripts import taxonomy
from nose.plugins.attrib import attr
from sequana import sequana_data


class TestPipeline(object):

    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.prog = "sequana_taxonomy"
        klass.params = {'prog': klass.prog}


    def test_help(self):
        try:
            taxonomy.main([self.prog, '--help', '1>/tmp/out', '2>/tmp/err'])
            assert False
        except SystemExit:
            pass
        else:
            raise Exception

    def test_input(self):
        try:
            df = taxonomy.main([self.prog, '--download', 'toydb'])
        except SystemExit:
            pass

    def test_analysis(self):
        file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
        file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz")

        # Test that database must be provided
        try:
            df = taxonomy.main([self.prog, '--file1', file1])
            assert False
        except:
            assert True

        from tempfile import TemporaryDirectory
        directory = TemporaryDirectory()
        df = taxonomy.main([self.prog, '--file1', file1, "--database", "toydb",
                "--file2", file2, "--verbose", "--output-directory", directory.name])
        from sequana import logger
        logger.info(directory.name)
        # cleanup

    def _test_wrong_db(self):
        try:
            df = taxonomy.main([self.prog, "--database", "dummy"])
            assert False
        except:
            assert True
