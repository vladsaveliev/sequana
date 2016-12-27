from sequana import sequana_data
from sequana.expdesign import ExpDesignAdapter
from sequana.expdesign import ExpDesignHiSeq, ExpDesignMiSeq


def test_designHiSeq():
    tt = ExpDesignHiSeq(sequana_data("test_expdesign_hiseq.csv", "testing"))
    assert list(tt.df['Index1_Seq'].values)  == ['TTAGGC', None, 'ACTTGA']


def test_designMiSeq():
    filename = sequana_data("test_expdesign_miseq_illumina.csv")
    tt = ExpDesignMiSeq(filename)
    tt.df.Index1_ID[0] == 1
    assert tt.adapter_type == "NEXTFlex-PCRfree"


def test_designMiSeq2():
    filename = sequana_data("test_expdesign_miseq_illumina2.csv")
    tt = ExpDesignMiSeq(filename)
    tt.df.Index2_Seq[0] == "ACGTCTCG"


def test_design_constructor():
    filename = sequana_data("test_expdesign_miseq_illumina.csv")
    tt = ExpDesignMiSeq(filename)

    # using existing design
    tt = ExpDesignAdapter(tt)
    # or a filename
    tt = ExpDesignAdapter(filename)

    # constructor for hiseq
    tt = ExpDesignAdapter(sequana_data("test_expdesign_hiseq.csv"))

    tt = ExpDesignAdapter(sequana_data("test_expdesign_generic.csv"))
    tt
    print(tt)
