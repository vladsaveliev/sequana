from sequana.bamtools import is_bam, is_sam, is_cram
from sequana.fastq import is_fastq
from sequana.fasta import is_fasta

def sniffer(filename):

    try:
        if is_sam(filename): return "SAM"
    except:
        pass

    try:
        if is_bam(filename): return "BAM"
    except:
        pass

    try:
        if is_cram(filename): return "CRAM"
    except:
        pass

    try:
        if is_fastq(filename): return "FASTQ"
    except:
        pass

    try:
        if is_fasta(filename): return "FASTA"
    except:
        pass

