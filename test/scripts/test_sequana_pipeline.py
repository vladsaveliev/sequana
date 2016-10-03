from sequana.scripts import main
from nose.plugins.attrib import attr
import mock
from sequana import sequana_data
import os


#@attr('skip')
class TestPipeline(object):

    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.prog = "sequana"
        klass.params = {'prog': klass.prog}

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""
        try:os.remove('quality.rules')
        except:pass
        try:os.remove('config.yaml')
        except:pass

        import shutil
        try:shutil.rmtree("Hm2_test")
        except:pass

        try:shutil.rmtree("report")
        except:pass

    def test_version(self):
        main.main([self.prog, '--version'])

    def test_init(self):
        try:
            # py3
            with mock.patch('builtins.input', return_value="y"):
                try:
                    main.main([self.prog, '--pipeline', "qualitydummy"])
                    assert False
                except:
                    assert True
        except:
            # py2
            with mock.patch('__builtin__.input', return_value="y"):
                main.main([self.prog, '--pipeline', "quality"])

            with mock.patch('__builtin__.input', return_value="y"):
                try:
                    main.main([self.prog, '--pipeline', "qualitydummy"])
                    assert False
                except:
                    assert True

    def test_run(self):
        pass

    def test_help(self):
        try:
            main.main([self.prog, '--help'])
            assert False
        except SystemExit:
            pass
        else:
            raise Exception

    @attr('onweb')
    def test_info(self):
        main.main([self.prog, '--info', "quality"])

    def test_show_pipelines(self):
        main.main([self.prog, '--show-pipelines'])

    def test_mutually_exclusive(self):
        try:
            main.main([self.prog, '--pipeline', 'quality', '--info'])
            assert False
        except:
            assert True

    def test_input(self):
        file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
        file2 = sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz', 'data')
        main.main([self.prog, "--pipeline", "quality", "--file1", file1, "--file2", file2, "--project", "Hm2_test"])




