from sequana import vcf_filter
from . import data
pathdata = data.__path__[0]
from easydev import TempFile
import os


def test_vcf_filter():
    ft = TempFile()
    vcf_output_expected = pathdata + os.sep + "vcf_filter_output.vcf"
    vcf_record = vcf_filter.VCF(pathdata + os.sep + "vcf_filter_input.vcf")
    filter_dict = {"QUAL": 1000, "FREQ": 90, "INFO": {"DP": ">10"}}
    vcf_record.filter_vcf(filter_dict, ft.name)
    flag = True
    with open(ft.name, 'r') as fl1:
        with open(vcf_output_expected, 'r') as fl2:
            for line1 in fl1:
                line2 = fl2.readline()
                if(line1 != line2):
                    flag = False
                    break
    assert flag == True
    ft.delete()
