import filecmp

from easydev import TempFile
from sequana import sequana_data
from sequana.freebayes_vcf_filter import VCF_freebayes


def test_vcf_filter():
    vcf_output_expected = sequana_data('JB409847.expected.vcf')
    v = VCF_freebayes(sequana_data('JB409847.vcf'))
    filter_dict = {'freebayes_score': 200, 'frequency': 0.85, 'min_depth': 10,
                   'forward_depth': 3, 'reverse_depth': 3, 'strand_ratio': 0.3}
    filter_v = v.filter_vcf(filter_dict)
    with TempFile(suffix='.vcf') as ft:
        filter_v.to_vcf(ft.name)
        compare_file = filecmp.cmp(ft.name, vcf_output_expected)
        assert compare_file

def test_constructor():
    try:
        VCF_freebayes('dummy')
        assert False
    except FileNotFoundError:
        assert True

def test_to_csv():
    filter_dict = {'freebayes_score': 200, 'frequency': 0.85, 'min_depth': 10,
                   'forward_depth': 3, 'reverse_depth': 3, 'strand_ratio': 0.3}
    v = VCF_freebayes(sequana_data('JB409847.expected.vcf'))
    filter_v = v.filter_vcf(filter_dict)
    with TempFile(suffix='.csv') as ft:
        filter_v.to_csv(ft.name)
