from sequoia import snakemake


def test_rules():
    "dag" in snakemake.rules.keys()
