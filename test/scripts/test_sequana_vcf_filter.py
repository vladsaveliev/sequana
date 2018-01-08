from sequana.scripts import vcf_filter
from sequana import sequana_data
from easydev import TempFile

prog = "sequana_vcf_filter"

def test_input():
    filename = sequana_data('test_vcf_mpileup_4dot1.vcf')

    with TempFile() as fout1:
        with TempFile() as fout2:

            res = vcf_filter.main([prog, '--input', filename, "--quality", "0",
                "--filter", "DP4[2]<2", "--filter",  "DP<50", 
                "--output", fout1.name,
                "--output-filtered", fout2.name])

            assert res == {'filtered': 209, 'N': 573, 'unfiltered': 364}

