from sequana import snaketools, sequana_data
from sequana.snaketools import DOTParser
import os
from nose.plugins.attrib import attr


testing = lambda x: sequana_data(x, "testing")

def test_dot_parser():
    s = DOTParser(testing("test_dag.dot"))
    s.add_urls()
    try:os.remove("test_dag.ann.dot")
    except:pass


def test_modules():
    assert "dag" in snaketools.modules.keys()
    assert snaketools.modules['dag'].endswith("dag.rules")


def test_getcleanup_rules():
    filename =  snaketools.modules['fastq_sampling']
    try:
        snaketools.get_cleanup_rules(filename)
    except:
        pass


def test_snakemake_stats():
    # this is created using snakemake with the option "--stats stats.txt"
    s = snaketools.SnakeMakeStats(testing("test_snakemake_stats.txt"))
    s.plot()


def test_module():
    # a rule without README
    m = snaketools.Module('mark_duplicates')
    m.description

    # a rule with README
    m = snaketools.Module('dag')
    m.description
    m.overview
    m.is_executable()
    m.check()

    # a pipeline
    m = snaketools.Module('quality_control')
    m.is_executable()
    m.check()
    m.snakefile
    m.name

@attr("onweb")
def test_module_onweb():
    m = snaketools.Module('quality_control')
    m.onweb()

def test_valid_config():
    config = snaketools.SequanaConfig(None)

    s = snaketools.Module("quality_control")
    config = snaketools.SequanaConfig(s.config)

    from easydev import TempFile
    with TempFile() as fh:
        config.save(fh.name)

def test_sequana_config():
    s = snaketools.Module("quality_control")
    config = snaketools.SequanaConfig(s.config)

    #assert config.config.get("kraken:kraken_database_directory") == '%(kraken_database_directory)s'
    #assert config.config.get("kraken:kraken_database_directory", "test") == "test"
    assert config.config.get("kraken:dummy", "test") == "test"
    assert config.config.get("kraken:dummy") == None


def test_file_name_factory():
    import glob

    def inner_test(ff):
        len(ff)
        print(ff)
        ff.filenames
        ff.realpaths
        ff.all_extensions
        ff.pathnames
        ff.extensions

    #list
    list_files = glob.glob("*.py")
    ff = snaketools.FileFactory(list_files)
    inner_test(ff)

    # glob
    ff = snaketools.FileFactory("*py")
    inner_test(ff)
    

    directory = os.path.dirname(sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz"))

    ff = snaketools.FastQFactory(directory + "/*fastq.gz", verbose=True)
    assert ff.tags == ['Hm2_GTGAAA_L005']

    ff.get_file1(ff.tags[0])
    ff.get_file2(ff.tags[0])
    assert len(ff) == 1








