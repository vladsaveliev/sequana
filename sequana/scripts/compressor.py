# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Standalone dedicated to taxonomic content (kraken)"""
import os
import sys
import shutil
import sys
from optparse import OptionParser
import argparse

from snakemake import shell
#from easydev import TempFile
import tempfile

from sequana.scripts.tools import SequanaOptions


class Options(argparse.ArgumentParser, SequanaOptions):
    def  __init__(self, prog="sequana_compressor"):
        usage = """Welcome to SEQUANA - Fastq compression standalone

    This standalone fetches recursively all files in a given format (--source)
    and transform them into another format (--target)

    Supported files must have one of the following extension:

        - fastq
        - fastq.gz
        - fastq.bz2

    The underlying compression tools used are pigz and pbzip2, which must be
    installed.

    sequana_compressor --source fastq.gz   --target fastq.bz2
    sequana_compressor --source fastq      --target fastq.bz2
    sequana_compressor --source fastq.gz   --target fastq
    sequana_compressor --source fastq.bz2  --target fastq


AUTHORS: Thomas Cokelaer
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana

        """
        description = """DESCRIPTION:
        """
        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("--source", dest="source", type=str,
            help="""Search for all source files with this extension. Valid
extensions are: fastq, fastq.gz, fastq.bz2, fastq.dscr, fq, fq.gz, fq.bz2 and
fq.dsrc""")
        self.add_argument("--target", dest="target", type=str,
            help="""Convert the source files to a new format. Valid
extensions are: fastq, fastq.gz, fastq.bz2, fastq.dscr, fq, fq.gz, fq.bz2 and
fq.dsrc""")
        self.add_argument("--recursive", dest="recursive",
            default=False,
            action="store_true", help="""recursive search""")
        self.add_argument("--threads", dest="threads",
            default=4,
            help="""number of threads to use per job (4)""")
        self.add_argument("--jobs", dest="jobs",
            default=4,
            help="""number of jobs to use at most (4) """)
        self.add_argument("--snakemake-options", dest="snakemake",
            default="--keep-going",
            help="""any valid list of options accepted by snakemake except
            -s and -j and -p""")
        self.add_version(self)
        self.add_verbose(self)
        self.add_cluster(self)

def main(args=None):

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
       options = user_options.parse_args(args[1:])

    if options.version:
        import sequana
        print(sequana.version)
        sys.exit()
    # valid codecs:
    valid_extensions = [("fastq." + ext2).rstrip(".")
                        for ext2 in ['', 'bz2', 'gz', 'dsrc']]

    valid_extensions += [("fq." + ext2).rstrip(".")
                        for ext2 in ['', 'bz2', 'gz', 'dsrc']]

    valid_combos = [(x, y) for x in valid_extensions
                           for y in valid_extensions
                           if x!=y]

    if (options.source, options.target) not in valid_combos:
        raise ValueError("""--target and --source combo not valid.
Must be in one of fastq, fastq.gz, fastq.bz2 or fastq.dsrc""")

    from easydev import TempFile
    # Create the config file locally
    with TempFile(suffix=".yaml", dir=".") as temp:
        fh = open(temp.name, "w")
        fh.write("compressor:\n")
        fh.write("    source: %s\n" %options.source)
        fh.write("    target: %s\n" % options.target)
        fh.write("    threads: %s\n" % options.threads)
        fh.write("    recursive: %s\n" % options.recursive)
        fh.write("    verbose: %s\n" % options.verbose)
        fh.close() # essential to close it because snakemake will try to use seek()
        from sequana import Module
        rule = Module("compressor").path + os.sep +  "compressor.rules"

        cmd = 'snakemake -s %s  --configfile %s -j %s -p ' % \
                (rule, temp.name, options.jobs)

        if options.cluster:
            cluster = ' --cluster "%s" ' % options.cluster
            cmd += cluster

        if options.snakemake:
            cmd += options.snakemake

        if options.verbose:
            print(cmd)

        shell(cmd)
        fh.close()


if __name__ == "__main__":
   import sys
   main(sys.argv)

