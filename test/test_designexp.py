from sequana import sequana_data
from sequana.designexp import ExpDesignHiSeq, ExpDesignMiSeq


def test_designHiSeq():
    tt = ExpDesignHiSeq(sequana_data("test_expdesign_hiseq.csv", "testing"))
    assert list(tt.df['Index Seq'].values)  == ['TTAGGC', 'GATCAG', 
        'ACTTGA', 'TAGCTT', 'ACTTGA']


def test_designMiSeq():
    tt = ExpDesignMiSeq(sequana_data("test_expdesign_miseq_illumina.csv", "testing"))
    tt.df.I7_Index_ID[0] == "NF01"
