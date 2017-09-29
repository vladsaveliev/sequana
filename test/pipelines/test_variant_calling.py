from sequana import SequanaConfig, sequana_data
from easydev import shellcmd
import subprocess
import json
import os
import tempfile
from .common import Pipeline

class VariantCallingPipeline(Pipeline):

    def __init__(self, wk=None):
        super(VariantCallingPipeline, self).__init__(wk=wk)
        # Define the data
        data = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
        input_directory = os.path.dirname(data)
        self.input_pattern = input_directory + "/Hm*gz"
        self.pipeline = "variant_calling"

        # Define the project and config file
        subprocess.check_call([
            "sequana", "--pipeline", self.pipeline,
            "--input-pattern", '%s' % self.input_pattern,
            "--working-directory", self.wk, "--force"
            ])
        # Add reference in the config
        cfg = SequanaConfig(self.wk + "/config.yaml")
        # We added a TTTT in position 5881
        cfg._yaml_code['bwa_mem_ref']['reference'] = sequana_data("measles.fa")
        cfg.save(self.wk + '/config.yaml')

    def check(self):
        from sequana.freebayes_vcf_filter import VCF_freebayes, Variant
        vcf = VCF_freebayes(self.wk +
            "/report_vc_Hm2_GTGAAA_L005/outputs/Hm2_GTGAAA_L005.raw.vcf")
        vcf.rewind()
        vv = [Variant(v)._resume for v in vcf]
        assert len(vv) == 85
        vv[29] == {'alternative': 'G',
            'chr': 'ENA|K01711|K01711.1',
            'depth': 5,
            'freebayes_score': 126.901,
            'frequency': '1.00',
            'position': '276',
            'reference': 'C',
            'strand_balance': '0.40'}


def test_variant_calling():
    QC = VariantCallingPipeline()
    try:
        QC.run()
        QC.check()
        QC.clean()
    except:
        QC.clean()

