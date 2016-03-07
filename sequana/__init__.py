__version__ = "$Rev: 10 $"
import pkg_resources
try:
    version = pkg_resources.require("sequana")[0].version
except:
    version = __version__





from .snakemake import rules
from .snakemake import ValidateConfig

from .samtools import BAM, SAM
