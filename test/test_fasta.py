import os

from sequana import FastA, sequana_data
from easydev import TempFile


def test_format_contigs_denovo():
    # test with a custom fasta
    filename = sequana_data("test_fasta.fasta") 
    contigs = FastA(filename)
    with TempFile(suffix='.fasta') as fh:
        contigs.format_contigs_denovo(fh.name)
