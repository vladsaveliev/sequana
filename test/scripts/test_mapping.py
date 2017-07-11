from sequana.scripts import mapping
from sequana import sequana_data
import os
import pytest
import shutil

prog = "sequana_mapping"


def test_analysis():
    file1 = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz")
    file2 = sequana_data("Hm2_GTGAAA_L005_R2_001.fastq.gz")
    reference = sequana_data("measles.fa")

    from tempfile import TemporaryDirectory
    directory = TemporaryDirectory()
    shutil.copy(file1, directory.name)
    shutil.copy(file2, directory.name)
    shutil.copy(reference, directory.name)

    df = mapping.main([prog, 
            '--file1', file1,
            "--file2", file2, 
            "--reference", reference])


def test_help():
    try:
        mapping.main([prog, '--help', '1>/tmp/out', '2>/tmp/err'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception
