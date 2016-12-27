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
""".. rubric:: Standalone application dedicated to coverage"""
import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse
from argparse import RawTextHelpFormatter


from sequana import bedtools
from sequana.reporting import report_mapping
from sequana.reporting import report_chromosome
from sequana.reporting import report_main

from easydev import shellcmd

from pylab import show, figure, savefig

from easydev.console import purple

# http://stackoverflow.com/questions/18462610/argumentparser-epilog-and-description-formatting-in-conjunction-with-argumentdef
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_coverage"):
        usage = purple("""\nWelcome to SEQUANA -- Coverage standalone

    Extract and plot coverage of one or more chromosomes/contigs in a BED or BAM
    file. In addition, running median used in conjunction with double thresholds
    extract regions of interests (low or high coverage). A reference may be
    provided to plot the coverage versus GC content. 

    The input file should be one of the following:

    - a BED file that is a tabulated file at least 3 columns.
      The first column being the reference, the second is the position 
      and the third column contains the coverage itself. 
    - or a BAM file that is converted automatically
      into a BED file using the following command:

        samtools depth -aa input.bam > output.bed

    If the reference is provided, an additional plot showing the coverage versus
    GC content is also shown.

    Here are some examples

        sequana_coverage --input file.bed --window-median 1001
        sequana_coverage --input file.bam --window-median 1001 -r <REFERENCE.fa>

    An other interesting option is to provide a BED file with 4 columns. The
    fourth column being another coverage data created with a filter. One can
    create such a file only from the BAM file using samtools as follows:

        samtools view -q 35  -o data.filtered.bam input.fa.sorted.bam
        samtools depth data.filtered.bam input.fa.sorted.bam  -aa > test.bed
        sequana_coverage --input test.bed --show-html 

    Note that the first file is the filtered one, and the second file is the
    unfiltered one.

        """)

        epilog = purple("""
----

AUTHORS: Thomas Cokelaer, Dimitri Desvillechabrol
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """)

        description = """DESCRIPTION:
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description, epilog=epilog,
                formatter_class=CustomFormatter)

        # options to fill the config file
        group = self.add_argument_group("Required argument")
        group.add_argument("-i", "--input", dest="input", type=str,
            help=("Input file in BED or BAM format. If a BAM file is "
                 "provided, it will be converted locally to a BED file "
                 "using genomecov, which must be installed."))

        group = self.add_argument_group("Optional biological arguments")
        group.add_argument('-c', "--chromosome", dest="chromosome", type=int,
            default=1,
            help=(  "Chromosome number (if only one, no need to use: the first"
                    " and only chromosome is chosen automatically). Default is"
                    " first chromosome found in the BED file. You may want to"
                    " analyse all chromosomes at the same time. If so, set this"
                    " parameter to -1"))
        group.add_argument('-o', "--circular", dest="circular",
            default=False, action="store_true",
            help="""If the DNA of the organism is circular (typically
            viruses or bacteria), set to True""")

        group = self.add_argument_group("General")
        group.add_argument("--output-directory", dest="output_directory",
            default="report", help="name of the output (report) directory.")
        group.add_argument('--show', dest="show", default=False,
            action='store_true', help="""Show the pictures (matplotlib)""")
        group.add_argument('--show-html', dest="show_html", default=False,
            action='store_true',
            help="""When report is created, you can open
            the main page automatically with this option (default is False)""")
        group.add_argument("-q", "--quiet", dest="verbose",
            default=True, action="store_false")
        group.add_argument('--no-report', dest="create_report",
            default=True, action='store_false',
            help="""Do not create any HTML report""")

        group = self.add_argument_group('Annotation')
        group.add_argument("-b", "--genbank", dest="genbank",
            type=str, default=None, help='a valida genbank annotation')

        # Group related to GC content
        group = self.add_argument_group("GC content related")
        group.add_argument('-r', "--reference", dest="reference", type=str,
            default=None,
            help="""If available, you can provide a reference (ENA/NCBI). It
                 must have the same length as the one used to create the
                 BAM or BED file. If provided, it is used to create the
                 coverage versus GC content image""")
        group.add_argument("-g", "--window-gc", dest="w_gc", type=int,
            help="""Length of the running window to compute the GC content""", default=200)
        group.add_argument('-n', "--nlevels", dest="levels", type=int,
            default=3, help="""Number of levels in the contour""")

        #group running median
        group = self.add_argument_group("Running Median related")
        group.add_argument("-w", "--window-median", dest="w_median", type=int,
            help="""Length of the running median window (default 4001,
                 recommended for viruses).  For long genome, 20001
                 or 30001 is recommended but larger windows may be
                 useful in the presence of long deleted regions.""",
            default=4001)

        group.add_argument("-k", "--mixture-models", dest="k", type=int,
            help="""Number of mixture models to use (default 2, although if sequencing
        depth is below 8, k is set to 1 automatically). To ignore that behaviour
        set k to the required value""",
            default=None)

        group.add_argument("-L", "--low-threshold", dest="low_threshold",
            default=None, type=float,
            help=("lower threshold (zscore) of the confidence interval. "
                "Overwrite value given by --threshold/-T"))
        group.add_argument("-H", "--high-threshold", dest="high_threshold",
            default=None, type=float,
            help=("higher threshold (zscore) of the confidence interval. "
                "Overwrite value given by --threshold/-T"))
        group.add_argument("-T", "--threshold", dest="threshold",
            default=4, type=float,
            help="""set lower and higher thresholds of the confidence interval.""")

        group = self.add_argument_group("Download reference")
        group.add_argument("--download-reference", dest="accession",
            default=None, type=str)
        group.add_argument("--download-genbank", dest="download_genbank",
            default=None, type=str)
        group.add_argument("--database", dest="database",
            default="ENA", type=str,
            choices=["ENA", "EUtils"],
            help="Download the reference from one of these database (default ENA)")


def main(args=None):

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
        options = user_options.parse_args(args[1:])

    if options.accession:
        print("Download accession %s from %s\n" %
            (options.accession, options.database))

        from bioservices.apps import download_fasta as df
        df.download_fasta(options.accession, method=options.database)
        return

    if options.download_genbank:
        from sequana.snpeff import download_fasta_and_genbank
        download_fasta_and_genbank(options.download_genbank,
                                   options.download_genbank, 
                                   genbank=True, fasta=False)
        return

    # We put the import here to make the --help faster
    from sequana import GenomeCov

    if options.genbank:
        assert os.path.exists(options.genbank), \
            "%s does not exists" % options.genbank

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

    if options.low_threshold is None:
        options.low_threshold = -options.threshold

    if options.high_threshold is None:
        options.high_threshold = options.threshold

    gc = GenomeCov(bedfile, options.low_threshold, options.high_threshold, 
            0.5, 0.5)

    if options.reference:
        print('Computing GC content')
        gc.compute_gc_content(options.reference, options.w_gc)

    if len(gc.chr_list) == 1:
        if options.verbose:
            print("There is only one chromosome. Selected automatically.")
        chrom = gc.chr_list[0]
        run_analysis(gc, chrom, 1, options) 
    elif options.chromosome <-1 or options.chromosome > len(gc.chr_list) or\
            options.chromosome == 0:
        raise ValueError("invalid --chromosome value ; must be in [1-%s]" % len(gc.chr_list)+1)
    else:
        # For uses, we start at position 1 but in python, we start at zero
        if options.chromosome == -1:
            chromosomes = [this for this in range(len(gc.chr_list)) ]
        else:
            chromosomes = [options.chromosome-1]


        N = len(chromosomes)
        for i,chrom_index in enumerate(chromosomes):
            chrom = gc.chr_list[chrom_index]
            chrom_name = gc.chr_list[chrom_index].chrom_name
            if options.verbose:
                print("There are %s chromosomes/contigs." % len(gc.chr_list))
                print("==================== analysing chrom/contig %s/%s (%s)" % (i+1,N,chrom_name))
            run_analysis(gc, chrom, chrom_index, options)


def run_analysis(gc, chrom, chrom_index, options):

    if options.verbose:
        print(chrom)

    if options.verbose:
        print('Computing running median')

    chrom.running_median(n=options.w_median, circular=options.circular)

    stats = chrom.get_stats(output="dataframe")
    stats.set_index("name", inplace=True)

    DOC = stats.ix['DOC'].Value
    if options.k is None and DOC < 8:
        options.k = 1
    elif options.k is None:
        options.k = 2

    if options.verbose:
        print("Number of mixture model %s " % options.k)
        print('Computing zscore')
    chrom.compute_zscore(k=options.k, verbose=options.verbose)

    if options.verbose:
        print("Computing centralness")

    # Let us save the thresholds first and then change it to compute centralness
    thresholds = chrom.thresholds.copy()

    chrom.thresholds.low = -3
    chrom.thresholds.high = 3
    c3 = chrom.get_centralness()

    chrom.thresholds.low = -4
    chrom.thresholds.high = 4
    c4 = chrom.get_centralness()
    chrom.thresholds = thresholds.copy()   # Get back to the original values

    if options.verbose and chrom.thresholds:
        print(chrom.thresholds)

    if options.verbose:
        res = chrom._get_best_gaussian()
        print("sigma and mu of the central distribution: mu=%s, sigma=%s" %
(round(res["mu"],3), round(res['sigma'],3)))
        print("Evenness: %8.3f" % chrom.get_evenness())
        print("Centralness (3 sigma): %f" % round(c3,3))
        print("Centralness (4 sigma): %f" % round(c4,4))

    if options.verbose:
        print("\n\n")

    figure(1)
    chrom.plot_coverage()

    # With the reference, we can plot others plots such as GC versus coverage
    if options.reference:
        figure(2)
        chrom.plot_gc_vs_coverage(Nlevels=options.levels, fontsize=20)
        filename_gc_cov = "coverage_vs_gc.chrom{0}.png".format(chrom_index)
        savefig(filename_gc_cov)

    if options.show:
        show()

    # Report chromosomes
    if options.verbose:
        print("Creating report in %s" % options.output_directory)


    if options.genbank:
        if options.verbose:
            print('Genbank: %s' % options.genbank)
        from sequana.tools import genbank_features_parser
        features = genbank_features_parser(options.genbank)
    else:
        features = None

    if options.create_report:
        sample_name = "coverage"
        r = report_chromosome.ChromosomeMappingReport(chrom,
                directory=options.output_directory, project="coverage", 
                sample=sample_name, features=features, verbose=options.verbose)

        if options.reference:
            from snakemake import shell
            shell("cp %s report/images" %  filename_gc_cov)
            r.jinja['coverage_vs_gc'] = """<img src="images/%s">""" % filename_gc_cov

        r.jinja['standalone_command'] = " ".join(sys.argv)
        r.create_report()
        if options.verbose:
            print("Report created. See ./report directory content and look for the HTML file. You can also use --show-html option")
        if options.show_html:
            from easydev import onweb
            onweb("report/%s_mapping.chrom%s.html" % (sample_name,
                  chrom_index))


if __name__ == "__main__":
   import sys
   main()#sys.argv)

