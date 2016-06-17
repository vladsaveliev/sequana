# Import -----------------------------------------------------------------------

from sequana import sequana_data
from sequana import snpeff
from easydev import TempFile

# Test -------------------------------------------------------------------------

def test_snpeff():

    # a custom refrence
    mydata = snpeff.SnpEff(reference=sequana_data("test_snpeff_ref.gb"))

    # an existing reference
    mydata = snpeff.SnpEff(
            reference="Bordetella_pertussis_Tohama_I_uid57617")

    with TempFile() as fh:
        mydata.launch_snpeff(sequana_data("test.vcf"), fh.name)

