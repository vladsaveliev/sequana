import os
import tempfile

from sequana.reporting.report_cutadapt import CutAdaptReport
from sequana import bedtools, sequana_data


def test_report():
    fh = tempfile.TemporaryDirectory()
    car = CutAdaptReport("sample_name", direcotry=fh.name)

    filename = sequana_data("test_cutadapt_paired.txt")
    car.read_data(filename)
    car.parse()
    car.get_histogram_data()
