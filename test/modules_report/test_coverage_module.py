from sequana import bedtools, sequana_data
from sequana.modules_report.coverage import CoverageModule
from sequana.utils import config


def test_coverage_module(tmpdir):

    bed = bedtools.GenomeCov(sequana_data("JB409847.bed"))
    fasta = sequana_data("JB409847.fasta")
    bed.compute_gc_content(fasta)
    c = bed.chr_list[0]
    c.run(4001)

    directory = tmpdir.mkdir('test_coverage_module')
    config.output_dir = str(directory)
    config.sample_name = "JB409847"
    CoverageModule(bed)
