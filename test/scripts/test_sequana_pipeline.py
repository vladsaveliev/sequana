from sequana.scripts import main
from sequana import sequana_data
import os


prog = "sequana"

"""    def teardown_class(klass):
        try:os.remove('quality.rules')
        except:pass
        try:os.remove('config.yaml')
        except:pass

        import shutil
        try:shutil.rmtree("Hm2_test")
        except:pass

        try:shutil.rmtree("report")
        except:pass
"""

def test_version():
    main.main([prog, '--version'])

def test_init(mocker):
    try:
        # py3
        with mocker.patch('builtins.input', return_value="y"):
            try:
                main.main([prog, '--pipeline', "qualitydummy"])
                assert False
            except:
                assert True
    except:
        # py2
        with mocker.patch('__builtin__.input', return_value="y"):
            main.main([prog, '--pipeline', "quality_control"])
        with mocker.patch('__builtin__.input', return_value="y"):
            try:
                main.main([prog, '--pipeline', "qualitydummy"])
                assert False
            except:
                assert True

def test_help():
    try:
        main.main([prog, '--help'])
        assert False
    except SystemExit:
        pass
    else:
        raise Exception

def test_info():
    main.main([prog, '--info', "quality"])

def test_show_pipelines():
    main.main([prog, '--show-pipelines'])

#def test_mutually_exclusive():
#    try:
#        main.main([prog, '--pipeline', 'quality', '--info'])
#        assert False
#    except:
#        assert True

def test_input():
    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    file2 = sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz', 'data')
    main.main([prog, "--pipeline", "quality_control", "--file1", 
              file1, "--file2", file2, "--no-adapters" , "--force"])




