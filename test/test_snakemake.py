from sequana import snakemake
import os
from . import data
pathdata = data.__path__[0]


def test_rules():
    "dag" in snakemake.rules.keys()



def test_snakemake_stats():

    # this is created using snakemake with the option "--stats stats.txt"
    s = snakemake.SnakeMakeStats(pathdata + os.sep + "stats.txt")
    s.plot()
    assert 'dag' in s.parse_data()['rules']

