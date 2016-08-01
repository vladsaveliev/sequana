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
"""Standalone application dedicated to coverage"""
import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse


from sequana import bedtools
from sequana.reporting import report_mapping
from sequana.reporting import report_chromosome
from sequana.reporting import report_main

from easydev import shellcmd

class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_coverage"):
        usage = """Welcome to SEQUANA - Coverage standalone

            sequana_coverage --input file.bed --window-median 1001
            sequana_coverage --input file.bam --window-median 1001 -r <REFERENCE.fa>

AUTHORS: Thomas Cokelaer, Dimitri Desvillechabrol
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        group = self.add_argument_group("Required argument")
        group.add_argument("-i", "--input", dest="input", type=str,
            required=True, help="Input file in BED or BAM format")

        group = self.add_argument_group("Optional biological arguments")
        group.add_argument('-c', "--chromosome", dest="chromosome", type=int,
            default=1,
            help="""Chromosome number (if only one, no need to use)""") 
        group.add_argument('-o', "--circular", dest="circular",
            default=False, action="store_true", help="""""") 

        group = self.add_argument_group("General")
        group.add_argument('--show', dest="show", default=False,
            action='store_true')
        group.add_argument('--show-html', dest="show_html", default=False,
            action='store_true')
        group.add_argument("-q", "--quiet", dest="verbose", 
            default=True, action="store_false")

        group = self.add_argument_group("GC content related")
        group.add_argument('-r', "--reference", dest="reference", type=str,
            default=None,help="""reference""") 
        group.add_argument("-g", "--window-gc", dest="w_gc", type=int,
            help="""Length of the running median window""", default=200)
        group.add_argument('-n', "--nlevels", dest="levels", type=int,
            default=3,
            help="""Number of levels in the contour""") 

        #group running median
        group = self.add_argument_group("Running Median related")
        group.add_argument("-w", "--window-median", dest="w_median", type=int,
            help="""Length of the running median window""", default=1001)
        
        group.add_argument("-k", "--mixture-models", dest="k", type=int,
            help="""Number of mixture models to use (default 2). If DOC is below
2, k is set to 1). To ignore that behavious set k to the required value""",
            default=None)

        group.add_argument("-L", "--low-threshold", dest="low_threshold",
            default=-3, type=float,
            help="lower threshold (zscore) of the confidence interval")
        group.add_argument("-H", "--high-threshold", dest="high_threshold",
            default=3, type=float,
            help="higher threshold (zscore) of the confidence interval")
        group.add_argument("-T", "--threshold", dest="threshold",
            default=3, type=float,
            help="higher threshold of the confidence interval")



def main(args=None):

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
        options = user_options.parse_args(args[1:])

    # We put the import here to make the --help faster
    from sequana import GenomeCov
    if options.verbose:
        print("Reading %s" % options.input)

    if options.input.endswith(".bam"):
        bedfile = options.input.replace(".bam", ".bed")
        if options.verbose:
            print("Converting BAM into BED file")
        shellcmd("bedtools genomecov -d -ibam %s > %s" % (options.input, bedfile))
    elif options.input.endswith(".bed"):
        bedfile = options.input
    else:
        raise ValueError("Input file must be a BAM or BED file")

    gc = GenomeCov(bedfile)

    if options.reference:
        print('Computing GC content')
        gc.compute_gc_content(options.reference, options.w_gc)

    if len(gc.chr_list) == 1:
        chrom = gc.chr_list[0]
    elif options.chromosome < 0 or options.chromosome > len(gc.chr_list):
        raise ValueError("invalid --chromosome value ; must be in [1-%s]" % len(gc.chr_list)+1)
    else:
        chrom = gc.chr_list[options.chromosome-1]

    if options.verbose:
        print('Computing running median')

    chrom.running_median(n=options.w_median, circular=options.circular)

    if options.verbose:
        print('Computing zscore')

    DOC = chrom.get_stats()['DOC']
    if options.k is None and DOC < 8:
        options.k = 1
    elif options.k is None:
        options.k = 2
    print("DOC: %s " % DOC)
    print("Number of mixture model %s " % options.k)


    chrom.compute_zscore(k=options.k)

    from pylab import show, figure, savefig

    figure(1)
    chrom.plot_coverage(low_threshold=options.low_threshold,
        high_threshold=options.high_threshold)

    # With the reference, we can plot others plots such as GC versus coverage
    if options.reference:
        figure(2)
        chrom.plot_gc_vs_coverage(Nlevels=options.levels, fontsize=20)
        savefig("coverage_vs_gc.png")

    if options.show:
        show()

    # Report mapping
    #r = report_mapping.MappingReport(directory="report",
    #    project="coverage")
    #r.set_data(gc)
    #r.create_report()

    if options.threshold:
        options.low_threshold = -options.threshold
        options.high_threshold = options.threshold

    # Report chromosomes
    chrom_index = options.chromosome
    print("Creating report")
    r = report_chromosome.ChromosomeMappingReport(chrom_index,
        low_threshold=options.low_threshold, 
        high_threshold=options.high_threshold,
        directory="report", project="coverage")
    r.set_data(chrom)
    r.create_report()
    if options.show_html:
        from easydev import onweb
        onweb("report/coverage_mapping.chrom%s.html" % chrom_index)


if __name__ == "__main__":
   import sys
   main()#sys.argv)

