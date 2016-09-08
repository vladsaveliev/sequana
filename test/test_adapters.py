from sequana import adapters
from sequana import sequana_data, FastA
from easydev import TempFile



def test_fasta_fwd_rev_to_columns():
    a1 = sequana_data("adapters_PCR-free_PF1_220616_fwd.fa", "data/adapters")
    a2 = sequana_data("adapters_PCR-free_PF1_220616_rev.fa", "data/adapters")
    f1 = FastA(a1)
    f2 = FastA(a2)
    assert f1 == f1
    assert f1 != f2
    assert len(f1) == 49
    assert len(f2) == 49

    with TempFile() as fh:
        adapters.fasta_fwd_rev_to_columns(a1, a2, fh.name)
    with TempFile() as fh:
        adapters.fasta_fwd_rev_to_columns(a1, None, output_filename=fh.name)
    with TempFile() as fh:
        adapters.fasta_fwd_rev_to_columns(a1, a2)
    with TempFile() as fh:
        adapters.fasta_fwd_rev_to_columns(a1, None)


def test_clean_ngs():
    a1 = sequana_data("adapters_PCR-free_PF1_220616_fwd.fa", "data/adapters")
    with TempFile() as fh:
        adapters.adapters_to_clean_ngs(a1, fh.name)


def test_adapter():
    a1 = adapters.Adapter("test", "ACGT", "comment")
    a2 = adapters.Adapter("test", "ACGT", "comment")
    assert a1 == a2


def test_adapters_removal_parser():
    data = sequana_data("test_adapter_removal_output.txt", "testing")
    results = adapters.adapter_removal_parser(data)
    assert sorted(results.keys()) == ["adapter1", "adapter2"]


def test_adapter_reader():
    from sequana.adapters import AdapterReader as AR
    data = sequana_data("adapters_with_duplicates.fa", "testing")
    try:
        AR(data)
    except ValueError:
        pass


    data1 = sequana_data("adapters_Nextera_PF1_220616_fwd.fa", "data/adapters")
    data2 = sequana_data("adapters_Nextera_PF1_220616_rev.fa", "data/adapters")
    data3 = sequana_data("adapters_Nextera_PF1_220616_revcomp.fa", "data/adapters")

    # try different constructors
    ar1 = AR(data1)
    ar_same = AR(ar1._data)    # from a list of dictionaries
    assert ar1 == ar_same
    ar_same = AR(ar1)           # from a AR instance
    assert ar1 == ar_same

    # __eq__
    assert len(ar1) == 56

    # accessors
    ar1.sequences, ar1.identifiers, ar1.comments

    ar1.get_adapter_by_sequence("ACGT")
    assert ar1.get_adapter_by_index("dummy") is None
    assert ar1.get_adapter_by_identifier("Nextera_index_N517|index_dna:N517")

    ar2 = AR(data2)
    ar2.reverse()
    # fails due to S516 ????????
    assert ar1 == ar2

    ar3 = AR(data3)
    ar3.reverse_complement()
    assert ar1 == ar3



def test_find_adapters_from_index_mapper():
    from sequana.adapters import FindAdaptersFromIndex
    ad = FindAdaptersFromIndex(sequana_data("test_index_mapper.csv", "testing"),
            "Nextera")
    assert ad.get_adapters("C4405-M1-EC1")
    ad.sample_names

    fwd, rev = ad.save_adapters_to_fasta("C4405-M1-EC1")
    import os
    os.remove(fwd)
    os.remove(rev)
