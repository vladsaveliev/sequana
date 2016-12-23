# coding: utf-8
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
import os
import shutil
import glob
import sys
import argparse

import pandas as pd


class Options(argparse.ArgumentParser):
    def __init__(self, prog="sequana_reports"):
        usage = """Welcome to SEQUANA - Reports generator
            
            sequana_reports --modules freebayes --file variants.vcf

AUTHORS: Thomas Cokelaer, Dimitri Desvillechabrol
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:

        Create HTML reports for differents results files.

        You can analyse these differents files:
            
            - csv from sequana_coverage
            - vcf from freebayes
        """
        super().__init__(usage=usage, prog=prog, description=description)

        # Options to create all requested report and generated


def main(args=None):
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
        options = user_options.parse_args(args[1:])
    options.verbose = not options.quiet

