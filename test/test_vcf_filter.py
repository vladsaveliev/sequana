import os
import filecmp

from easydev import TempFile
from sequana import vcf_filter
from . import data
pathdata = data.__path__[0]


def test_vcf_filter():
    ft = TempFile()
    vcf_output_expected = pathdata + os.sep + "vcf_filter_output.vcf"
    vcf_record = vcf_filter.VCF(pathdata + os.sep + "vcf_filter_input.vcf")
    filter_dict = {"freebayes_score": 10000, "frequency": 0.85,
                   "min_depth": 10, "forward_depth": 3, "reverse_depth": 3,
                   "strand_ratio": 0.2}
    vcf_record.filter_vcf(filter_dict, ft.name)
    compare_file = filecmp.cmp(ft.name, vcf_output_expected)
    assert compare_file
    ft.delete()


def test_constructor():
    vcf_filter.VCF('dummy')
