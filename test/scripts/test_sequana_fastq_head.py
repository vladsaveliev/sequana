from sequana.scripts import fastq_head
from sequana import sequana_data


prog = "sequana_fastq_head"



def test_input():

    from easydev import TempFile
    fh = TempFile(suffix=".fastq.gz")
    filename = sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz')
    df = fastq_head.main([prog, '--input', filename, '--nlines', 
        "100", "--output", fh.name])

    df = fastq_head.main([prog, filename, "100", fh.name])
    fh.delete()




