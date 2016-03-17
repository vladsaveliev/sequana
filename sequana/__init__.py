__version__ = "$Rev: 10 $"
import pkg_resources
try:
    version = pkg_resources.require("sequana")[0].version
except:
    version = __version__





from .snakemake import rules
from .snakemake import ValidateConfig
from .phred import Quality
from .samtools import BAM, SAM
from .fastq import FastQ, FastQC




def sequana_data(filename, where=None):
    """Simple utilities to retrieve data sets from gdsctools/share directory"""
    import os
    import easydev
    share = easydev.get_shared_directory_path('sequana') 
    share = os.sep.join([share, 'data'])
    # in the code one may use / or \ 
    if where:
        filename = os.sep.join([share, where, filename])
    else:
        filename = os.sep.join([share, filename])
    if os.path.exists(filename) is False:
        raise Exception('unknown file %s' % filename)
    return filename

