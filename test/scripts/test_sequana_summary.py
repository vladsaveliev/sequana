from sequana.scripts import summary
from sequana import sequana_data



prog = "sequana_summary"

def _test_version():
    summary.main([prog, '--version'])

def test_help():
    try:
        summary.main([prog, '--help'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception

def test_input():
    filename = sequana_data('virus.bed', 'data')
    df = summary.main([prog, '--file', filename])
    len(df)

    filename = sequana_data('test.fastq', "testing")
    df = summary.main([prog, '--file', filename])




