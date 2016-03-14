from sequana import vcf_filter
from . import data
pathdata = data.__path__[0]
from easydev import TempFile
import os
import filecmp


def test_vcf_filter():
    ft = TempFile()
    vcf_output_expected = pathdata + os.sep + "vcf_filter_output.vcf"
    vcf_record = vcf_filter.VCF(pathdata + os.sep + "vcf_filter_input.vcf")
    filter_dict = {"QUAL": 10000, "FREQ": 0.85, "INFO": {"DP": ">10",
                                                        "AO": ">200",
                                                        "SRP": "<100"}}
    vcf_record.filter_vcf(filter_dict, ft.name)
    compare_file = filecmp.cmp(ft.name, vcf_output_expected)
    assert compare_file == True
    ft.delete()
