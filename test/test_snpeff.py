from sequana import sequana_data
from sequana import snpeff
from easydev import TempFile
import os

def test_snpeff():
    # a custom refrence
    mydata = snpeff.SnpEff(reference=sequana_data("JB409847.gbk"))
    with TempFile() as fh:
        mydata.launch_snpeff(sequana_data("JB409847.vcf"), fh.name) 

    # cleanup
    try:
        os.remove("snpEff.config")
    except:
        pass

    try:
        os.remove("snpEff_genes.txt")
    except:
        pass

    try:
        os.remove("snpEff_summary.html")
    except:
        pass



def test_snpeff_download():
    with TempFile() as fh:
        snpeff.download_fasta_and_genbank("K01711", fh.name)

    with TempFile() as fh:
        try:
            snpeff.download_fasta_and_genbank("dummyK01711", fh.name)
            assert False
        except ValueError:
            assert True
        except:
            assert False
