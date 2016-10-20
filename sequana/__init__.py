__version__ = "$Rev: 10 $"
import pkg_resources
try:
    version = pkg_resources.require("sequana")[0].version
except:
    version = __version__

from easydev import CustomConfig
configuration = CustomConfig("sequana", verbose=False)
sequana_config_path = configuration.user_config_dir

# This must be import before all other modules (sequana_data function)
from .datatools import sequana_data

# snakemake related
from .snaketools import modules
from .snaketools import SequanaConfig
from .snaketools import Module, ModuleFinder, FastQFactory

# various utilities for IO/data analysis
from .adapters import AdapterReader, FindAdaptersFromIndex, Adapter
from .bamtools import BAM, SAMFlags
from .bedtools import GenomeCov
from .coverage import Coverage
from .fastq import FastQ, FastQC, Identifier
from .fasta import FastA
from .kraken_builder import KrakenBuilder
from .krona import KronaMerger
from .kraken import KrakenResults, KrakenPipeline, KrakenAnalysis, KrakenDownload
from .phred import Quality
from .running_median import RunningMedian
from .snpeff import SnpEff
from .vcf_filter import VCF

# Reports
from sequana.reporting.report_bam import BAMReport
from sequana.reporting.report_fastqc import FastQCReport
from sequana.reporting.report_fastq_stats import FastQStatsReport
from sequana.reporting.report_summary import SequanaSummary

# The standalone app
from . import scripts


