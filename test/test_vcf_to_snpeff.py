# Import -----------------------------------------------------------------------

import os
import shutil
from sequana import vcf_to_snpeff
from easydev import TempFile
from . import data
pathdata = data.__path__[0]

# Test -------------------------------------------------------------------------

def test_vcf_to_snpeff():

    mydata = vcf_to_snpeff.Vcf_to_snpeff(
            vcf_filename=pathdata + os.sep + "test.vcf",
            reference=pathdata + os.sep + "test.gb")
    mydata.add_custom_db()
    mydata = vcf_to_snpeff.Vcf_to_snpeff(
            vcf_filename=pathdata + os.sep + "test.vcf",
            reference="Bordetella_pertussis_Tohama_I_uid57617")
    with TempFile(suffix='.png') as fh:
        mydata.launch_snpEff(fh.name)
