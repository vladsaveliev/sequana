from sequana.tools import bam_to_mapped_unmapped_fastq
from sequana import sequana_data


def test_bam2fastq():
    data = sequana_data("test.bam", "testing")
    res = bam_to_mapped_unmapped_fastq(data)
