import os
import tempfile
from sequana import sequana_data, SequanaConfig
import subprocess

from .common import Pipeline


class RNASeqPipeline(Pipeline):
    def __init__(self, wk=None):
        super(RNASeqPipeline, self).__init__(wk)
        data = sequana_data("KO_ATCACG_R1_test.fastq.gz")
        input_directory = os.path.dirname(data)
        self.input_pattern = input_directory + "/KO_ATCACG_R1_test.fastq.gz"
        self.pipeline = "rnaseq"

        #self.output = self.wk + "/Hm2_GTGAAA_L005/report_qc_Hm2_GTGAAA_L005/summary.json"
        subprocess.check_call([
            "sequana", "--pipeline", self.pipeline,
            "--input-pattern", '%s'% self.input_pattern,
            "--working-directory", self.wk,
            "--adapter-fwd", 
            "GATCGGAAGAGCACACGTCTGAACTCCAGTCA", 
            "--adapter-rev", 
            "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC",
            "--force"])

        # Need to edit the config file
        cfg = SequanaConfig(self.wk + "/config.yaml")
        cfg._yaml_code['genome']['genome_directory'] = "Saccer3"
        cfg._yaml_code['genome']['name'] = "Saccer3"
        cfg._yaml_code['genome']['fasta_file'] = "Saccer3/Saccer3.fa"
        cfg._yaml_code['genome']['fasta_file'] = "Saccer3/Saccer3.gff"
        cfg.save(self.wk + '/config.yaml')

        #TODO        Download the Saccer3 data

    def run(self):
        pass

    def check(self):
        pass


def test_pipeline():
    QC = RNASeqPipeline()
    try:
        QC.run()
        QC.check()
        QC.clean()
    except:
        QC.clean()
        raise Exception



