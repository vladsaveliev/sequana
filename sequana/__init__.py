__version__ = "$Rev: 10 $"
import pkg_resources
try:
    version = pkg_resources.require("sequana")[0].version
except:
    version = __version__


# snakemake related
from .snaketools import modules
from .snaketools import SequanaConfig
from .snaketools import Module, ModuleFinder


# tools
from .bamtools import BAM, SAMFlags
from .bedtools import Genomecov
from .coverage import Coverage
from .fastq import FastQ, FastQC, Identifier
from .fasta import FastA
from .kraken import KrakenResults, KrakenTaxon, KronaMerger, KrakenBuilder
from .phred import Quality
from .vcf_filter import VCF

# Reports
from .report_bam import BAMReport

# The standalone app
from . import scripts


def sequana_data(filename, where=None):
    """Simple utilities to retrieve data sets from gdsctools/share directory"""
    import os
    import easydev
    sequana_path = easydev.get_package_location('sequana')
    share = os.sep.join([sequana_path , "sequana", 'resources'])
    # in the code one may use / or \ 
    if where:
        filename = os.sep.join([share, where, filename])
    else:
        filename = os.sep.join([share, filename])
    if os.path.exists(filename) is False:
        raise Exception('unknown file %s' % filename)
    return filename

