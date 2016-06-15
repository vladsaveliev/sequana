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
from .kraken_builder import KrakenBuilder
from .krona import KronaMerger
from .kraken import KrakenResults, KrakenPipeline, KrakenAnalysis
from .phred import Quality
from .snpeff import SnpEff
from .vcf_filter import VCF

# Reports
from .report_bam import BAMReport

# The standalone app
from . import scripts


def sequana_data(filename=None, where=None):
    """Simple utilities to retrieve data sets from gdsctools/share directory"""
    import os
    import easydev
    import glob
    sequana_path = easydev.get_package_location('sequana')
    sharedir = os.sep.join([sequana_path , "sequana", 'resources'])


    if filename is None:
        for thisdir in ['data', 'testing']:
            for filename in glob.glob(sharedir + "/%s/*" % thisdir):
                to_ignore = ["__init__.py", "__pycache__"]
                if filename.endswith('.pyc') or filename in to_ignore:
                    pass
                else:
                    print('sequana("%s", "%s")' % (os.path.split(filename)[1], thisdir))
        raise ValueError("Choose a valid file from the list above")
    # in the code one may use / or \ 
    if where:
        filename = os.sep.join([sharedir, where, filename])
    else:
        filename = os.sep.join([sharedir, filename])

    if os.path.exists(filename) is False:
        raise Exception('unknown file %s. Type sequana_data() to get a list of valid names' % filename)


    return filename

