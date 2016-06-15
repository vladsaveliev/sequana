# Import -----------------------------------------------------------------------

import os
import shutil
from sequana import snpeff
from easydev import TempFile
from . import data
pathdata = data.__path__[0]

# Test -------------------------------------------------------------------------

def test_snpeff():

    mydata = snpeff.SnpEff(reference=pathdata + os.sep + "test.gb")
    mydata = snpeff.SnpEff(
            reference="Bordetella_pertussis_Tohama_I_uid57617")
    with TempFile() as fh:
        mydata.launch_snpeff(pathdata + os.sep + "test.vcf", fh.name)
