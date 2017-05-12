from sequana.sequence import DNA, RNA, Repeats
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


def test_repeats():
    rep = Repeats(datafile)
    rep.threshold = 11
    assert len(rep.begin_end_repeat_position) == 158
    rep.do_merge = True
    assert len(rep.begin_end_repeat_position) == 156
    rep.do_merge = False
    assert len(rep.begin_end_repeat_position) == 158
    assert rep.df_shustring.shape == (15888,2)
    assert rep.longest_shustring == 15

    # test histogram
    rep = Repeats(datafile)
    rep.threshold = 11
    rep.hist_length_repeats()



def test_gc_skew():

    data = sequana_data("measles.fa")
    dna = DNA(data)
    dna.window = 100
    dna.plot_all_skews()
