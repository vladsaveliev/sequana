import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse



class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_coverage"):
        usage = """Welcome to SEQUANA - Coverage standalone

            sequana_coverage --bed file.bed --window-size 1001

AUTHORS: Thomas Cokelaer, Dimitri Desvillechabrol
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("-b", "--bed", dest="bedfile", type=str,
            required=True, help="""filename of a BED file""")
        self.add_argument("-w", "--window-median", dest="w_median", type=int,
            help="""Length of the running median window""", default=1001)
        self.add_argument("-g", "--window-gc", dest="w_gc", type=int,
            help="""Length of the running median window""", default=200)
        self.add_argument('-c', "--chromosome", dest="chromosome", type=int,
            default=1,
            help="""Chromosome number (if only one, no need to use)""") 
        self.add_argument('-n', "--nlevels", dest="levels", type=int,
            default=3,
            help="""Number of levels in the contour""") 
        self.add_argument('-r', "--reference", dest="reference", type=str,
            default=None,help="""reference""") 

        self.add_argument('-o', "--circular", dest="circular",
            default=False, action="store_true", help="""""") 


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
    print("Reading %s" % options.bedfile)
    gc = GenomeCov(options.bedfile)
    if options.reference:
        print('Computing GC content')
        gc.compute_gc_content(options.reference, options.w_gc)


    print(len(gc.chr_list))
    if len(gc.chr_list) == 1:
        gc = gc.chr_list[0]
    elif options.chromosome < 0 or options.chromosome > len(gc.chr_list):
        raise ValueError("invalid --chromosome value ; must be in [1-%s]" % len(gc.chr_list)+1)
    else:
        gc = gc.chr_list[options.chromosome-1]


    print('Computing running median')
    gc.running_median(n=options.w_median, circular=options.circular)
    print('Computing zscore')
    gc.compute_zscore()

    from pylab import show, figure

    figure(1)
    gc.plot_coverage()
    if options.reference:
        figure(2)
        gc.plot_gc_vs_coverage(Nlevels=options.levels, fontsize=20)
    show()


if __name__ == "__main__":
   import sys
   main()#sys.argv)

