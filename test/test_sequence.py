from sequana.sequence import DNA, RNA
from sequana import sequana_data


datafile = sequana_data("measles.fa")


def test_dna():
    seq = "ACGTTTT"
    dna = DNA(seq)
    assert dna.sequence == seq
    assert dna.get_complement() == "TGCAAAA"
    assert dna.get_reverse() == "TTTTGCA"
    assert dna.get_reverse_complement() == "AAAACGT"
    dna.check()
    dna.stats()

    # inplace functions
    dna.reverse()
    dna.complement()
    dna.reverse_complement()
 

    dna = DNA("jjjj")
    try:
        dna.check()
        assert False
    except:
        assert True


    # read a file and tests the __len__ method
    dna = DNA(datafile)
    assert len(dna) == 15894
    dna.gc_content()

    # test occurences
    dna._data = "ACGTGGGGGTT"
    assert dna.get_occurences("GGG", False) == [4]
    assert dna.get_occurences("GGG", True) == [4, 5, 6]


def test_rna():
    rna = RNA("ACUG")
