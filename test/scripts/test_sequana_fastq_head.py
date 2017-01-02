from sequana.scripts import fastq_head
from nose.plugins.attrib import attr
from sequana import sequana_data


#@attr("skip")
class TestPipeline(object):

    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        klass.prog = "sequana_fastq_head"
        klass.params = {'prog': klass.prog}



    def test_input(self):

        from easydev import TempFile
        fh = TempFile(suffix=".fastq.gz")
        filename = sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz')
        df = fastq_head.main([self.prog, '--input', filename, '--nlines', 
            "100", "--output", fh.name])

        df = fastq_head.main([self.prog, filename, "100", "--output"])
        fh.delete()




