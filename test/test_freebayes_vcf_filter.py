import filecmp

from easydev import TempFile
from sequana import sequana_data
from sequana.freebayes_vcf_filter import VCF_freebayes


def test_vcf_filter():
    vcf_output_expected = sequana_data('JB409847.expected.vcf')
    vcf_record = VCF_freebayes(sequana_data('JB409847.vcf'))
    filter_dict = {'freebayes_score': 200, 'frequency': 0.85, 'min_depth': 10,
                   'forward_depth': 3, 'reverse_depth': 3, 'strand_ratio': 0.3}
    with TempFile(suffix='.vcf') as ft:
        vcf_record.filter_vcf(filter_dict, ft.name)
        compare_file = filecmp.cmp(ft.name, vcf_output_expected)
        assert compare_file

def test_constructor():
    VCF_freebayes('dummy')

def test_to_csv():
    filter_dict = {'freebayes_score': 200, 'frequency': 0.85, 'min_depth': 10,
                   'forward_depth': 3, 'reverse_depth': 3, 'strand_ratio': 0.3}
    vcf_record = VCF_freebayes(sequana_data('JB409847.expected.vcf'))
    with TempFile(suffix='.csv') as ft:
        vcf_record.to_csv(ft.name, filter_dict)
