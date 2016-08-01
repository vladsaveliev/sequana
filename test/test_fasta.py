import os

from sequana import FastA, sequana_data


def test_format_contigs_denovo():
    # test with a custom fasta
    filename = sequana_data("test_fasta.fasta") 
    contigs = FastA(filename)
    contigs.format_contigs_denovo(project="test")

    # cleanup
    try:
        os.remove("test.ab500.fasta")
        os.remove("test.bl500.fasta")
    except:
        pass
