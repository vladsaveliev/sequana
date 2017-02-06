from sequana.scripts import fastq_count
from sequana import sequana_data


prog = "sequana_fastq_count"

def test_input():
    filename = sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz')
    df = fastq_count.main([prog, '--input', filename])




