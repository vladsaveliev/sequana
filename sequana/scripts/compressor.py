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

import tempfile

from sequana.scripts.tools import SequanaOptions
from sequana import Module, SequanaConfig


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


    If your job(s) were interrupted (ctrl+C), your directories will are
    problably locked. Use the --unlock option in such situations.

        sequana_compressor --source ... --target ... --unlock

    --source and --target must be provided but no analysis will be performed 
    at that stage.

    Then, type your command again.

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
            default=4, type=int,
            help="""Maximum number of threads to use per task (4).""")
        self.add_argument("--jobs", dest="jobs",
            default=4,
            help="""Maximum number of cores to use at most (4). """)

        self.add_argument("--unlock", action="store_true",
            help="""If you stopped the application, the underlying snakemake
                process are interrupted and directories were snakemake was
                launch will be locked. If so please use this option using the
                --source and --target as when you got the error message""")

        self.add_argument("--snakemake-options", dest="snakemake",
            default=" --keep-going ",
            help="""any valid list of options accepted by snakemake except
            -s and -j . Note that by default --keep-going is used ; If you set
            this argument yourself, you have to add --keep-going as well otherwise it stops
            at the first error encountered""")
        self.add_version(self)
        self.add_cluster(self)
        self.add_quiet(self)


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
Must be one of fastq, fastq.gz, fastq.bz2 or fastq.dsrc""")

    from easydev import TempFile
    # Create the config file locally

    module = Module("compressor")

    with TempFile(suffix=".yaml", dir=".") as temp:
        cfg = SequanaConfig(module.config)
        cfg.config.compressor.source = options.source
        cfg.config.compressor.target = options.target
        cfg.config.compressor.recursive = options.recursive
        cfg.config.compressor.verbose = options.verbose
        cfg.config.compressor.threads = options.threads
        cfg._update_yaml()
        cfg.save(filename=temp.name)


        # The Snakefile can stay in its original place:
        rule = module.path + os.sep +  "compressor.rules"

        # Run the snakemake command itself.
        cmd = 'snakemake -s %s  --configfile %s -j %s ' % \
                (rule, temp.name, options.jobs)

        if options.verbose is False:
            cmd += " --quiet "
        else:
            cmd += " -p "

        # for slurm only: --cores-per-socket
        if options.cluster:
            cluster = ' --cluster "%s" ' % options.cluster
            cmd += cluster

        if options.snakemake:
            cmd += options.snakemake

        if options.unlock:
            cmd += " --unlock "

        if options.verbose:
            print(cmd)

        # On travis, snakemake.shell command from snakemake fails.
        # Most probably because travis itself uses a subprocess.
        # excute from easydev uses pexpect.spawn, which seems to work well
        from easydev import execute
        execute(cmd, showcmd=False)


if __name__ == "__main__":
   import sys
   main(sys.argv)

