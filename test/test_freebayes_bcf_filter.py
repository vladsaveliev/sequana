import filecmp

from easydev import TempFile
from sequana import sequana_data
from sequana.freebayes_bcf_filter import BCF_freebayes


def test_bcf_filter():
    vcf_output_expected = sequana_data('JB409847.filter.vcf')
    bcf = BCF_freebayes(sequana_data('JB409847.bcf'))
    filter_dict = {'freebayes_score': 200, 'frequency': 0.85, 'min_depth': 10,
                   'forward_depth': 3, 'reverse_depth': 3, 'strand_ratio': 0.3}
    filter_bcf = bcf.filter_bcf(filter_dict)
    with TempFile(suffix='.vcf') as fp:
        filter_bcf.to_vcf(fp.name)
        compare_file = filecmp.cmp(fp.name, vcf_output_expected)
        assert compare_file

def test_constructor():
    try:
        BCF_freebayes('dummy')
        assert False
    except OSError:
        assert True 

def test_to_csv():
    filter_dict = {'freebayes_score': 200, 'frequency': 0.85, 'min_depth': 10,
                   'forward_depth': 3, 'reverse_depth': 3, 'strand_ratio': 0.3}
    bcf = BCF_freebayes(sequana_data('JB409847.bcf'))
    filter_bcf = bcf.filter_bcf(filter_dict)
    with TempFile(suffix='.csv') as ft:
        filter_bcf.to_csv(ft.name)
