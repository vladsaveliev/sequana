from sequana.tools import bam_to_mapped_unmapped_fastq, reverse_complement, StatsBAM2Mapped
from sequana import sequana_data
from sequana.tools import bam_get_paired_distance


def test_StatsBAM2Mapped():
    data = sequana_data("test.bam", "testing")
    res = StatsBAM2Mapped(data)


def test_bam2fastq():
    data = sequana_data("test.bam", "testing")
    res = bam_to_mapped_unmapped_fastq(data)



def test_reverse_complement():
    assert reverse_complement("AACCGGTTA") == 'TAACCGGTT'
def test_reverse():
    from sequana.tools import reverse
    assert reverse("AACCGG") == 'GGCCAA'


def test_distance():
    data = sequana_data("test.bam", "testing")
    distances = bam_get_paired_distance(data)


def test_gc_content():
    from sequana.tools import gc_content
    data = sequana_data('test.fasta', "testing")
    gc_content(data, 10)['seq1']
