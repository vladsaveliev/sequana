from sequana import adapters 
from sequana import sequana_data, FastA
from easydev import TempFile



def test_fasta_fwd_rev_to_columns():
    a1 = sequana_data("adapters_netflex_pcr_free_1_fwd.fa", "data")
    a2 = sequana_data("adapters_netflex_pcr_free_1_rev.fa", "data")
    f1 = FastA(a1)
    f2 = FastA(a2)
    assert f1 == f1
    assert f1 != f2
    assert len(f1) == 50
    assert len(f2) == 50

    with TempFile() as fh:
        adapters.fasta_fwd_rev_to_columns(a1, a2, fh.name)
    with TempFile() as fh:
        adapters.fasta_fwd_rev_to_columns(a1, None, output_filename=fh.name)
    with TempFile() as fh:
        adapters.fasta_fwd_rev_to_columns(a1, a2)
    with TempFile() as fh:
        adapters.fasta_fwd_rev_to_columns(a1, None)



def test_clean_ngs():
    a1 = sequana_data("adapters_netflex_pcr_free_1_fwd.fa", "data")
    with TempFile() as fh:
        adapters.adapters_to_clean_ngs(a1, fh.name)


def test_adapters_removal_parser():
    data = sequana_data("test_adapter_removal_output.txt", "testing")
    results = adapters.adapter_removal_parser(data)
    assert sorted(results.keys()) == ["adapter1", "adapter2"]


def test_adapters_db():

    a1 = sequana_data("adapters_netflex_pcr_free_1_fwd.fa", "data")
    a2 = sequana_data("adapters_netflex_pcr_free_1_rev.fa", "data")
    db = adapters.AdapterDB(a1)
    db.load_fasta(a2)
    assert len(db.df) == 100
    assert db.get_name(100000070) == "NextFlex_PCR_Free_adapter20_r"


    db = adapters.AdapterDB()
    db.load_all()


