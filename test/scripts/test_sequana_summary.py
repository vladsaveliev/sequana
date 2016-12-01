from sequana.scripts import summary
from nose.plugins.attrib import attr
from sequana import sequana_data


#@attr("skip")
class TestPipeline(object):

    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.prog = "sequana_summary"
        klass.params = {'prog': klass.prog}

    def _test_version(self):
        summary.main([self.prog, '--version'])

    def test_help(self):
        try:
            summary.main([self.prog, '--help'])
            assert False
        except SystemExit:
            pass
        else:
            raise Exception

    def test_input(self):
        filename = sequana_data('virus.bed', 'data')
        df = summary.main([self.prog, '--file', filename])
        len(df)

        filename = sequana_data('test.fastq', "testing")
        df = summary.main([self.prog, '--file', filename])




