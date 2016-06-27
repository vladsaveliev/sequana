from sequana.scripts import main
from nose.plugins.attrib import attr
import mock


class TestPipeline(object):

    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.prog = "sequana"
        klass.params = {'prog': klass.prog}

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""
        import os
        # local nosetests execution
        try:os.remove('README')
        except:pass
        try:os.remove('quality.rules')
        except:pass
        try:os.remove('config.yaml')
        except:pass


    def setUp(self):
        """This method is run once before _each_ test method is executed"""

    def teardown(self):
        """This method is run once after _each_ test method is executed"""

    def test_version(self):
        main.main([self.prog, '--version'])

    def test_init(self):
        try:
            # py3
            with mock.patch('builtins.input', return_value="y"):
                main.main([self.prog, '--pipeline', "quality", "--force-init"])

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

        #main.main([self.prog, '--pipeline', "quality", "--force-init"])

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


