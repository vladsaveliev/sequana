# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Extract head of a zipped or unzipped FastQ file"""
from sequana.vcf_filter import VCF

import sys
from optparse import OptionParser
import argparse
from sequana import logger
from easydev.console import purple

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


epilog = purple("""
----

AUTHORS: Thomas Cokelaer
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """)


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_vcf_filter"):
        usage = """%s Only for VCF using mpileup version 4.1 for now\n""" % prog
        usage += """usage2: %s vcf_filter""" % prog
        usage += """Examples:

    sequana_vcf_filter --input test.vcf --quality 40
                --filter "AF1>0.95&AF1<0.05"
                --filter "MQ<30"


    You should look into the VCF file to figure out the valid TAGs. Then, you
    can apply various filters. 

    A filter should be interpreted as :

    ''filter out variants that agree with the filter''

    For example::

        --filter "DP<30"

    means ''remove SNPs with DP below 30''. We accept those types of comparison:

        DP<30
        DP<=30
        DP>30
        DP>=30

    For some tags, you want to keep values within or outside a range of
    values. You can then use the & and | characters::

        DP<30|>60  # to keep only values in the ranges [0-30] and [60-infinite]

    or

        DP>30&<60  # to keep only values in the range [30-60]

    Some tags stores a list of values. For instance DP4 contains 4 values.
    To filter the value at position 1, use e.g.::

        DP4[0]<0.5

    you can use the same convention for the range as above::

        DP4[0]>0.05&<0.95

    you may also need something like:

        sum(DP4[2]+DP4[3]) <2


    Note that you must use quotes to surround the filter values.

    Instead of providing the filters one by one and forget about it, you can store
    them in a file with the following format (same for the --quality argument):

    [general]
    quality=50

    [filters]
    DP                  = <4
    DP4[0]              = <10
    DP4[2]              = <10
    sum(DP4[0], DP4[1]) = <10
    AF1                 = >0.05&<0.95

        """
        super(Options, self).__init__(usage=usage, prog=prog,
                epilog=epilog,
                formatter_class=CustomFormatter)

        self.add_argument("--input", dest='input_filename', type=str,
                            required=True, help="input fastq gzipped or not")

        self.add_argument("--quality", dest="quality",
                          type=int, default=0, help="filter sites below this quality")
        self.add_argument("--depth", dest="depth",
                          type=int, default=0, help="filter sites with depth below this number")

        self.add_argument("--filter", dest="filter", action="append",
                        nargs=1, type=str, default=[],
                        help="Provide as many filters as you want. See example above ")

        self.add_argument("--output", dest="output_filename",
                            default="remaining.vcf", type=str)

        self.add_argument("--output-filtered", dest="output_filtered_filename",
                            default="filtered.vcf", type=str)

        self.add_argument('--level', dest="level",
            default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

        self.add_argument("--filter-file", dest="filter_file", default=None,
                            help=r"""You may provide config file storing all
filters and quality (to re-use). See format explanation above.

""")


def main(args=None):
    if args is None:
        args = sys.argv[:]

    print("Welcome to sequana_vcf_filter")
    user_options = Options(prog="sequana_vcf_filter")

    if len(args) == 1 or "--help" in args:
        user_options.parse_args(["prog", "--help"])
    elif len(args) == 2:
        class SimpleOpt():
            pass
        options = SimpleOpt()
        options.input_filename = args[1]
    else:
        options = user_options.parse_args(args[1:])

    # set the level
    logger.level = options.level

    vcf = VCF(options.input_filename)
    vcf.vcf.filter_dict['QUAL'] =  options.quality
    vcf.vcf.filter_dict['INFO'] = {}

    print("------------------")
    # Read filters from a file
    if options.filter_file:
        import configparser
        cfg = configparser.RawConfigParser()
        cfg.optionsxform = str
        cfg.read(options.filter_file)
        if cfg.has_section('filters'):
            for key, value in cfg.items('filters'):
                vcf.vcf.filter_dict["INFO"][key.upper()] = value
        else:
            raise ValueError("filter file must contain a section "
                             "[filters] use --help for more information")

        if cfg.has_section('general'):
            quality = cfg.getint('general', 'quality')
            vcf.vcf.filter_dict['QUAL'] =  quality

    if options.quality != 0:
        vcf.vcf.filter_dict['QUAL'] =  options.quality


    for this in options.filter:
        this = this[0]
        signs = [">", "<", ">=", "<="]
        for sign in signs:
            if sign in this:
                key, value = this.split(sign, 1)
                key = key.strip()
                value = sign.strip() + value.strip()
                vcf.vcf.filter_dict['INFO'][key] = value
                break


    print(vcf.vcf.filter_dict)

    res = vcf.vcf.filter_vcf(options.output_filename,
                       output_filtered=options.output_filtered_filename)

    print()
    #print(res)
    return res


if __name__ == "__main__":
   import sys
   main(sys.argv)

