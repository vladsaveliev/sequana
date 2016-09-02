from sequana.sequence import DNA



def test_dna():
    seq = "ACGTTTT"
    dna = DNA(seq)
    assert dna.sequence == seq
    assert dna.get_complement() == "TGCAAAA"
    assert dna.get_reverse() == "TTTTGCA"
    assert dna.get_reverse_complement() == "AAAACGT"
    dna.check()
   

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
