from sequana import bedtools, sequana_data
from sequana.modules_report.coverage import CoverageModule
from sequana.utils import config


def test_coverage_module(tmpdir):
    bed = bedtools.GenomeCov(sequana_data("JB409847.cov.csv"),
                             sequana_data("JB409847.gbk"))
    directory = tmpdir.mkdir('test_coverage_module')
    config.output_dir = str(directory)
    config.sample_name = "JB409847"
    CoverageModule(bed)
