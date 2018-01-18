from sequana import Cigar



def test_cigar():

    cigar = "10S80M5I5I1D"  # the last 1D must be ignored, 5I5I should be understood
    c = Cigar(cigar)
    assert len(c) == 100

    assert c.as_tuple() == (("S",10),("M",80),("I",5), ("I",5), ("D",1))
    assert c.as_sequence() == "S"*10+"M"*80+"I"*10+"D"
    assert c.as_dict() == {"S":10, "M":80, "I":10, "D":1}

    print(c)
    c

    c = Cigar("1S1S10M")
    c.compress()
    assert c.cigarstring == "2S10M"

    c = Cigar("1S10M")
    c.compress()
    assert c.cigarstring == "1S10M"

    c = Cigar("1S1S1S1S")
    c.compress()
    assert c.cigarstring == "4S"
