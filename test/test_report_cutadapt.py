import os
import tempfile

from sequana.modules_report.cutadapt import CutadaptModule
from sequana import bedtools, sequana_data


def test_report():
    # Used the PCRFree adapters
    fh = tempfile.TemporaryDirectory()
    filename = sequana_data("test_cutadapt_paired.txt")
    mod = CutadaptModule(filename, "sample_name")
    assert mod.jinja['command'].startswith("cutadapt")
    assert mod.jinja['mode'] == 'Paired-end'
    assert mod.jinja['paired_reads1_with_adapters'] == '273'
    assert mod.jinja['paired_reads2_with_adapters'] == '243'
    assert mod.jinja['paired_reads_kept'] == '2,189'

    # atropos results should be identical except for a few differences (e.g.
    # command starts with atropos)
    # Note however that here, we used Nextera (I believe)
    filename = sequana_data("test_atropos_paired.txt")
    mod = CutadaptModule(filename, "sample_name")
    assert mod.jinja['command'].startswith("atropos")
    assert mod.jinja['mode'] == 'Paired-end'
    assert mod.jinja['paired_reads1_with_adapters'] == '197'
    assert mod.jinja['paired_reads2_with_adapters'] == '262'
    assert mod.jinja['paired_reads_kept'] == '2,420'
