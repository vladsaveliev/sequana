from sequana import sequana_data
from sequana.designexp import ExpDesignHiSeq, ExpDesignMiSeq


def test_designHiSeq():
    tt = ExpDesignHiSeq(sequana_data("test_expdesign_hiseq.csv", "testing"))
    assert list(tt.df['Index1_Seq'].values)  == ['TTAGGC', None, 'ACTTGA']


def test_designMiSeq():
    tt = ExpDesignMiSeq(sequana_data("test_expdesign_miseq_illumina.csv", "testing"))
    tt.df.Index1_ID[0] == 1


def test_designMiSeq2():
    tt = ExpDesignMiSeq(sequana_data("test_expdesign_miseq_illumina2.csv", "testing"))
    tt.df.Index2_Seq[0] == "ACGTCTCG"

