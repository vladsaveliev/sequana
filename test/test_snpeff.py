from sequana import sequana_data
from sequana import snpeff
from easydev import TempFile


def test_snpeff():

    # a custom refrence
    mydata = snpeff.SnpEff(reference=sequana_data("test_snpeff_ref.gb"))
    with TempFile() as fh:
        mydata.launch_snpeff(sequana_data("test.vcf"), fh.name)

