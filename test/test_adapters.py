import os
from sequana import adapters
from sequana import sequana_data, FastA
from easydev import TempFile
from sequana.adapters import FindAdaptersFromDesign


def test_fasta_fwd_rev_to_columns():
    a1 = sequana_data("adapters_PCRFree_fwd.fa")
    a2 = sequana_data("adapters_PCRFree_rev.fa")
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
    a1 = sequana_data("adapters_PCRFree_fwd.fa")
    with TempFile() as fh:
        adapters.adapters_to_clean_ngs(a1, fh.name)


def test_adapter():
    a1 = adapters.Adapter("test", "ACGT", "comment")
    a2 = adapters.Adapter("test", "ACGT", "comment")
    a3 = adapters.Adapter("test", "ACGC", "comment")
    a4 = adapters.Adapter("test", "ACGT", "comment2")
    assert a1 == a2
    assert a1 != a3
    assert a1 != a4

    try:
        # > should not be included
        adapters.Adapter(">test", "ACGC", "comment")
        assert False
    except ValueError:
        assert True

    a1 = adapters.Adapter("test", "ACGT", "comment")
    a1.sequence = "other"
    assert a1.sequence == "other"
    a1.identifier = "new"
    a1.comment = "new"

    # constructor with dictionary
    adapters.Adapter({"identifier":"test", "comment":"ACGT", "sequence":"ACFT"})

    # constructor with an adapter
    adapters.Adapter(a1)

    # check some attributes
    ad = adapters.Adapter("Nextera|name:N505|seq:ACGT","AAAACGTTTTT","comment")
    assert ad.name == 'N505'
    assert ad.index_sequence == 'ACGT'
    # case where there is no extra field
    ad = adapters.Adapter("Nextera|seq:ACGT","AAAACGTTTTT","comment")
    assert ad.name is None


    a1 = adapters.Adapter("id1", "ACGT", "test")
    a2 = adapters.Adapter("id2", "ACGT", "test")
    assert a1 != a2


def test_adapters_removal_parser():
    data = sequana_data("test_adapter_removal_output.txt", "testing")
    results = adapters.adapter_removal_parser(data)
    assert sorted(results.keys()) == ["adapter1", "adapter2"]


def test_adapter_reader():
    from sequana.adapters import AdapterReader as AR
    data = sequana_data("adapters_with_duplicates.fa", "testing")
    try:
        ar = AR(data)
        ar.sanity_check()
    except ValueError:
        pass

    data1 = sequana_data("adapters_Nextera_fwd.fa")
    data2 = sequana_data("adapters_Nextera_rev.fa")
    data3 = sequana_data("adapters_Nextera_revcomp.fa")

    # try different constructors
    ar1 = AR(data1)
    ar1.__repr__()
    ar1.index_sequences
    ar_same = AR(ar1._data)    # from a list of dictionaries
    assert ar1 == ar_same
    ar_same = AR(ar1)           # from a AR instance
    assert ar1 == ar_same
    assert ar1[0]['identifier'] == 'Universal_Adapter|name:universal'
    ar1.index_names
    assert ar1.get_adapter_by_sequence("XXX") is None
    try:
        ar1.get_adapter_by_identifier("XXX")
        assert False
    except ValueError:
        assert True

    # __eq__
    assert len(ar1) == 56

    # accessors
    ar1.sequences, ar1.identifiers, ar1.comments

    ar1.get_adapter_by_sequence("ACGT")
    assert ar1.get_adapter_by_index_name("dummy") is None
    assert ar1.get_adapter_by_identifier("Nextera_index_N517")

    ar2 = AR(data2)
    ar2.reverse()
    # fails due to S516 ????????
    assert ar1 == ar2

    ar3 = AR(data3)
    ar3.reverse_complement()
    assert ar1 == ar3

    # test to_fasta method
    with TempFile() as fh:
        ar1.to_fasta(fh.name)


def test_nextera():
    # simple indexing
    design = sequana_data("test_index_mapper.csv")
    ad = FindAdaptersFromDesign(design, "Nextera")
    results = ad.get_adapters_from_sample("C4405-M1-EC1")
    assert results['index1']['fwd'][0].identifier  == 'Nextera_index_N703|name:N703|seq:AGGCAGAA'

    ad.check()  # all samples are used in get_adapters_from_sample
    ad.sample_names

    fwd, rev = ad.save_adapters_to_fasta("C4405-M1-EC1")
    os.remove(fwd)
    os.remove(rev)

    # double indexing
    design = sequana_data("test_expdesign_hiseq_doubleindex.csv")
    fa = FindAdaptersFromDesign(design, "Nextera")


def test_pcrfree():
    design = sequana_data("test_index_mapper.csv")

    try:
        ad = FindAdaptersFromDesign(design, "error")
        assert False
    except Exception:
        assert True

    # Other input from PCRFree
    ad = FindAdaptersFromDesign(design, "PCRFree")

    # Test the index1/2_seq with 3 cases
    # index1 present only,
    # no index at all (None)
    # index1 and index2 present
    design1 = sequana_data("test_expdesign_hiseq.csv")
    ad1 = FindAdaptersFromDesign(design1, "PCRFree")
    ad1.check()
    res1 = ad1.get_adapters_from_sample("553-iH2-1")
    res2 = ad1.get_adapters_from_sample("539-st2")
    res3 = ad1.get_adapters_from_sample("107-st2")

    assert res1['index1']['fwd'].identifier == "NextFlex_PCR_Free_adapter8|name:8|seq:TTAGGC"
    assert res1['index1']['fwd'].name == "8"
    assert res1['index1']['rev'].name == "8"

    assert list(res2.keys()) == ["universal"]

    assert res3['index1']['fwd'].name == "9"
    assert res3['index1']['rev'].name == "9"
    assert res3['index2']['fwd'].name == "10"
    assert res3['index2']['rev'].name == "10"

    # double indexing
    # This is a double indexing for PCRFree, which has not been tested
    # since it requires 16S adapters not yet in sequana
    """design2 = sequana_data("test_expdesign_miseq_illumina2.csv")
    ad2 = FindAdaptersFromDesign(design2, "PCRFree")
    assert ad2.get_adapters_from_sample('M2')['index1']['fwd'].identifier == \
            'NextFlex_PCR_Free_adapter2|name:2|seq:TGACCA'
    assert ad2.get_adapters_from_sample('M2')['index2']['fwd'].identifier == \
            'NextFlex_PCR_Free_adapter13|name:13|seq:AGTCAA'
    """


    design = sequana_data("test_expdesign_miseq_illumina.csv")
    ad = FindAdaptersFromDesign(design, "PCRFree")
    res = ad.get_adapters_from_sample("CR81-L1236-P1")
    assert res['index1']['fwd'].identifier == 'NextFlex_PCR_Free_adapter1|name:1|seq:CGATGT'



    design1 = sequana_data("test_expdesign_miseq_illumina_1.csv")
    ad = FindAdaptersFromDesign(design1, "PCRFree")
    ad.check() # all sample names must be found
    res = ad.get_adapters_from_sample("CM-2685")['index1']['fwd']
    assert res.name == "3"


def test_wrong_design():
    design = sequana_data("test_expdesign_wrong.csv")
    ad = FindAdaptersFromDesign(design, "PCRFree")
    try:
        ad.check()
        assert False
    except:
        assert True


def test_rubicon():
    design = sequana_data("test_expdesign_rubicon.csv")
    fa = FindAdaptersFromDesign(design, "Rubicon")
    fa.check()


def test_get_sequana_adapters():

    assert "adapters_PCRFree_rev.fa" in adapters.get_sequana_adapters("PCRFree", "rev")
    assert "adapters_Rubicon_fwd.fa" in adapters.get_sequana_adapters("Rubicon", "fwd")
    assert "adapters_Nextera_revcomp.fa" in adapters.get_sequana_adapters("Nextera", "revcomp")

    try:
        adapters.get_sequana_adapters("Nexter", "fwd")
        assert False
    except ValueError:
        assert True

    try:
        adapters.get_sequana_adapters("Nextera", "fw")
        assert False
    except ValueError:
        assert True


# This also check the Small adapters
def test_duplicated_design():

    filename = sequana_data("test_expdesign_hiseq_duplicated_index.csv")
    ss = FindAdaptersFromDesign(filename, "Small")
    res = ss.get_adapters_from_sample("VB-22")
    assert res['index1']['fwd'].identifier == "Small_Adapter_5|name:small5|seq:ACAGTG"
    assert res['index1']['fwd'].sequence == "CAAGCAGAAGACGGCATACGAGATACAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
