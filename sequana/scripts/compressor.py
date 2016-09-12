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
import glob
import sys
from optparse import OptionParser
import argparse

from snakemake import shell
from easydev import TempFile 

class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_taxonomy"):
        usage = """Welcome to SEQUANA - Compression standalone

AUTHORS: Thomas Cokelaer
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("--from", dest="_from", type=str,
            help="""fastq, fastq.gz, fastq.bz2""")
        self.add_argument("--to", dest="_to", type=str,
            help="""fastq, fastq.gz, fastq.bz2 """)


def main(args=None):

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
       options = user_options.parse_args(args[1:])


    options.recursive = True

    # valid codecs:
    valid_combos = [
        ("fastq", "fastq.gz"),
        ("fastq", "fastq.bz2"),
        ("fastq.gz", "fastq"),
        ("fastq.bz2", "fastq"),
        ("fastq.bz2", "fastq.gz"),
        ("fastq.gz", "fastq.bz2")]

    # Create the config file
    temp = TempFile(suffix=".yaml")
    fh = open(temp.name, "w")
    fh.write("compressor:\n")
    fh.write("    source: %s\n" %options._from)
    fh.write("    target: %s\n" % options._to)
    fh.write("    recursive: %s\n" % options.recursive)
    fh.close() # essential to close it because snakemake will try to use seek()

    rule = "/home/cokelaer/Work/github/sequana/sequana/rules/compressor/"
    rule += "compressor.rules"

    cmd = "snakemake -s %s  --configfile %s -j 4 -p" % (rule, temp.name)
    shell(cmd)

    temp.delete()



if __name__ == "__main__":
   import sys
   main(sys.argv)

