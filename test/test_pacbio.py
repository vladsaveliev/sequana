from sequana.pacbio import BAMPacbio
from sequana import sequana_data

def test_pacbio():
    BAMPacbio(sequana_data("test_pacbio_subreads.bam"))
