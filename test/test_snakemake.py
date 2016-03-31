from sequana import snakemake
import os


def test_rules():
    "dag" in snakemake.rules.keys()



def test_snakemake_stats():

    # this is created using snakemake with the option "--stats stats.txt"
    s = snakemake.SnakeMakeStats(sequana_data("test_snakemake_stats.txt"))
    s.plot()
    assert 'dag' in s.parse_data()['rules']


def test_modules():
    m = snakemake.Modules()
    m.onweb('dag')
