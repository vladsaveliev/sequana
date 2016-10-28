from sequana.reporting.report_fastq_stats import FastQStatsReport
from sequana import sequana_data

filename = sequana_data("test_summary_fastq_stats.json")

def test_report():
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        report = FastQStatsReport(tmpdir)
        report.filenames = [filename]
        report.create_report()
