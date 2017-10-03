from sequana import SequanaConfig, sequana_data
from easydev import shellcmd
import subprocess
import json
import os
import tempfile
from .common import Pipeline


class PacbioQCPipeline(Pipeline):

    def __init__(self, wk=None):
        super(PacbioQCPipeline, self).__init__(wk=wk)
        # Define the data
        data = sequana_data("test_pacbio_subreads.bam")
        input_directory = os.path.dirname(data)
        self.input_pattern = input_directory + "/test_pacbio_subreads.bam"
        self.pipeline = "pacbio_qc"

        # Define the project and config file
        subprocess.check_call([
            "sequana", "--pipeline", self.pipeline,
            "--input-pattern", '%s' % self.input_pattern,
            "--extension", "bam",
            "--working-directory", self.wk, "--force",
            "--jobs", "1"
            ])

        cfg = SequanaConfig(self.wk + "/config.yaml")
        cfg._yaml_code["input_directory"] = ''
        cfg._yaml_code["input_readtag"] = "_R[12]_"
        cfg._yaml_code['input_extension'] = "bam"
        cfg._yaml_code['input_pattern'] = self.input_pattern 
        cfg._yaml_code['input_samples'] = "CommentedMap([('file1', None), ('file2', None)])"
        cfg.save(self.wk + '/config.yaml')

        self.output = self.wk + "/test_pacbio_subreads//summary_test_pacbio_subreads.json"


    def check(self):
        if os.path.exists(self.output):
            data = json.load(open(self.output))
        assert "hist_gc" in data




def test_pipeline():
    QC = PacbioQCPipeline()
    try:
        QC.run()
        QC.check()
        QC.clean()
    except:
        QC.clean()
        raise Exception
