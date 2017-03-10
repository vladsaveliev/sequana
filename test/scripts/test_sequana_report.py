import os

from sequana.scripts import reports
from sequana import sequana_data
from sequana.utils import config


prog = 'sequana_report'

def test_version():
    try:
        reports.main([prog, '--version'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception

def test_help():
    try:
        reports.main([prog, '--help'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception

def test_sequana_report(tmpdir):
    directory = tmpdir.mkdir('report')
    try:
        reports.main([prog, '--input-files', sequana_data('JB409847.cov.csv'),
                     '--output-directory', str(directory)])
        assert False
    except Exception as err:
        print(err)
        assert True
    assert os.path.exists(str(directory) + os.sep + 'sequana_coverage.html')        
