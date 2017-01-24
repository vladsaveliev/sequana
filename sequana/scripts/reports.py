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
import re
import json

import pandas as pd

from sequana.utils import config


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
        group = self.add_argument_group("Required arguments")
        group.add_argument(
            '-d', '--input-directory', type=str, dest='input_dir', help=(
            "Search for CSV and JSON files in the directory. It creates a "
            "report directory with a summary.html as main page."))
        group.add_argument('-i', '--input', type=str, nargs='*',
                           dest='input_file',
                           help="All CSV and JSON files produce by sequana")
        group.add_argument('-o', '--output-directory', type=str,
                           dest='output_dir', default='report',
                           help="Name of the output (report) directory")
        group = self.add_argument_group("General options")
        group.add_argument('--version', dest='version', action='store_true',
                           help="Print version of sequana")
        group.add_argument('-v', '--verbose', dest='verbose',
                           action='store_true', help="Display logs") 


def main(args=None):
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog='sequana')

    # If --help or no options provided, show the help
    if len(args) == 1:
        user_options.parse_args(['prog', '--help'])
    else:
        options = user_options.parse_args(args[1:])
    config.output_dir = options.output_dir

    # list of all target files (csv or json)
    input_files = [f for p in ('*.csv', '*.json') for f in
                   glob.glob(os.sep.join([options.input_dir, p]),
                   recursive=True)]

    for f in input_files:
        # get tools to retrieve module
        with open(f, 'r') as fp:
            if f.endswith('json'):
                data = json.load(fp)
                tool = data['tool']
            else:
                line = fp.readline()
                tool = re.findall("sequana_\w+", line)[0]
                data = f
        # load and execute module
        mod = config.module_dict[tool].load()
        mod(data)


if __name__ == '__main__':
    main()
