__version__ = "$Rev: 10 $"
import pkg_resources
try:
    version = pkg_resources.require("sequoia")[0].version
except:
    version = __version__




# a global function to get the templates from sequoia/share/templates
try:
    import os
    def get_template_path(name):
        # Is it a local directory ?
        if os.path.exists(name):
            return name
        else:
            import easydev
            template_path = easydev.get_shared_directory_path("sequoia")
            template_path += os.sep + "templates"  + os.sep + name
            return template_path
except:
    # before installation, this will fail but once installed,
    # template_path will be in hte namespace
    def get_template_path(name):
        return name



from .snakemake import rules
from .snakemake import ValidateConfig

from .samtools import BAM, SAM
