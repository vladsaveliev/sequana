from sequana.blaster import *
from sequana import sequana_data



def _test_blast_xml_parser():

    ebi = sequana_data('testing/test_blast_ebi.xml')
    local = sequana_data('testing/test_blast_local.xml')

    # read local files
    assert len(BlastHits(ebi).hits) == 389
    assert len(BlastHits(local).hits) == 6

    # read string
    assert len(BlastHits(open(ebi).read()).hits) == 389
    assert len(BlastHits(open(local).read()).hits) == 6


def _test_blast_online():

    bro = BlastRunOnline()
    results = bro.run("ACGTACGT", database="em_rel_vrl")
    # results should be empty but not guaranteed

    fastq = FastQ(sequana_data("testing/test.fastq.gz"))
    read = next(fastq)
    sequence = read['sequence']
    results = bro.run(sequence, database="em_rel_vrl")
    assert len(results.hits) > 0
