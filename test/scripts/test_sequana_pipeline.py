from sequana.scripts import main
from sequana import sequana_data
import pytest
from unittest.mock import patch

prog = "sequana"


def test_main():
    main.main([prog])


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

def test_info(mocker):
    def func(*args, **kwargs):
        pass
    with patch('easydev.onweb', func ):
        main.main([prog, '--info', "quality_control"])

# This mock does not work somehow, it opens the web page... is it in conflict
# with the test_info ?
def _test_issue(mocker):
    def func(*args, **kwargs):
        pass
    with patch('easydev.onweb', func ):
        main.main([prog, '--issue'])

def test_show_pipelines():
    main.main([prog, '--show-pipelines'])

def test_mutually_exclusive():
    try:
        main.main([prog, '--pipeline', 'quality_control', '--info', "quality_control"])
        assert False
    except SystemExit:
        assert True

def test_input(tmpdir):
    directory = tmpdir.mkdir("analysis")
    name = directory.__str__()

    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    file2 = sequana_data('Hm2_GTGAAA_L005_R2_001.fastq.gz', 'data')
    main.main([prog, "--pipeline", "quality_control", "--file1", 
              file1, "--file2", file2, "--no-adapters" , "--force",
              "--working-directory", name])

    try:
        main.main([prog, "--pipeline", "quality_control", "--file1", 
              file1+"dummy", "--no-adapters" , "--force",
              "--working-directory", name])
        assert False
    except ValueError:
        assert True
    try:
        main.main([prog, "--pipeline", "quality_control", "--file1", 
              file1, "--no-adapters" , "--force", "--file2",
              file2+"duly","--working-directory", name])
        assert False
    except ValueError:
        assert True
    try:
        main.main([prog, "--pipeline", "quality_control", "--file1", 
              file1, "--no-adapters" , "--force", "--file2",
              file2,"--working-directory", name, "--kraken", "dummy"])
        assert False
    except ValueError:
        assert True
    try:
        main.main([prog, "--pipeline", "quality_control", "--file1", 
              file1, "--no-adapters" , "--force", "--file2",
              file2,"--working-directory", name, "--kraken", "dummy"])
        assert False
    except ValueError:
        assert True


def test_without_cluster_config(tmpdir):
    directory = tmpdir.mkdir("analysis")
    name = directory.__str__()
    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    main.main([prog, "--pipeline", "rnaseq", "--file1", file1, "--snakemake-cluster", 
        '"sbatch --mem={cluster.ram}"', "--ignore-cluster-config", "--force",
        "--working-directory", name, "--no-adapters"])
    with open("%s/runme.sh" %  name) as fh: 
        assert "cluster_config" not in fh.read()

# Fails on travis
def _test_with_cluster_config(tmpdir):
    directory = tmpdir.mkdir("analysis2")
    name = directory.__str__()
    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    main.main([prog, "--pipeline", "rnaseq", "--file1", file1, "--snakemake-cluster", 
        '"sbatch --mem={cluster.ram}"',"--force",
        "--working-directory", name])
    with open("%s/runme.sh" %  name) as fh: 
        assert "cluster_config" in fh.read()


def test_get_config():

    main.main([prog,  "--pipeline", "quality_control",
        "--get-config"])
    import os
    os.remove("config.yaml")


def test_config_params():
    file1 = sequana_data('Hm2_GTGAAA_L005_R1_001.fastq.gz', 'data')
    main.main([prog, "--pipeline", "quality_control", "--force",  "--file1", file1,
        "--no-adapters",  "--config-params",  "bwa_mem_phix:threads:4"])

    try:
        main.main([prog, "--pipeline", "quality_control", "--force",  "--file1", file1,
            "--no-adapters",  "--config-params",  "bwa_mem_phix"])
        assert False
    except:
        assert True
