from sequana.modules_report.fastq_stats import FastQStatsModule
from sequana import sequana_data
import shutil
from sequana.utils import config

filename = sequana_data("test_summary_fastq_stats.json")

def test_report(tmpdir):
    directory = tmpdir.mkdir('test_module')
    shutil.copy(filename, str(directory))
    config.output_dir = str(directory) 

    #
    #report = FastQStatsModule(tmpdir, "dummy_link_fastqc")

    #report.create_report()
