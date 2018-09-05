from sequana.bamtools import is_bam, is_sam, is_cram


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


