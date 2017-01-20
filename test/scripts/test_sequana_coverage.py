from sequana.scripts import coverage
from sequana import sequana_data
import pytest

prog = "sequana_coverage"

@pytest.fixture
def coveragefix():
    import os
    # local nosetests execution
    try:os.remove('README')
    except:pass
    try:os.remove('quality.rules')
    except:pass
    try:os.remove('config.yaml')
    except:pass


def test_version():
    try:
        coverage.main([prog, '--version'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception


def test_help():
    try:
        coverage.main([prog, '--help'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception


def test_input(tmpdir):

    import os
    directory = tmpdir.mkdir("report")
    name = directory.__str__()

    filename = sequana_data('virus.bed', 'data')
    reference = sequana_data('tofill.fa', 'data')
    coverage.main([prog, '-i', filename, "-o", "--output-directory", name])
    assert os.path.exists(name + os.sep + "coverage_mapping.chrom1.html")
