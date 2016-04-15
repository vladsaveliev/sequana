from sequana.report_phix import PhixReport


def test_report():
    from sequana import sequana_data
    stats = sequana_data("test_snakemake_stats.txt", "testing")
    from sequana import Module
    module = Module('phix_removal')
    r = PhixReport()
    r.create_report()
