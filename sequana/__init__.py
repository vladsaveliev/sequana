
import pkg_resources
try:
    version = pkg_resources.require("sequana")[0].version
except:
    version = ">=0.18.0"

import colorlog as logger

from easydev import CustomConfig
configuration = CustomConfig("sequana", verbose=False)
sequana_config_path = configuration.user_config_dir

# This must be import before all other modules (sequana_data function)
from .datatools import sequana_data

from .snaketools import *
from .adapters import AdapterReader, FindAdaptersFromDesign, Adapter
from .expdesign import ExpDesignAdapter
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


# The standalone app
from . import scripts

