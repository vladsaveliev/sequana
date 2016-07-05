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
        self.add_argument("-w", "--window-size", dest="ws", type=int,
            help="""Length of the running median window""", default=1001)
        self.add_argument('-c', "--chromosome", dest="chromosome", type=int,
            default=1,
            help="""Chromosome number (if only one, no need to use)""") 


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
    from sequana import Genomecov
    print(options.bedfile)
    gc = Genomecov(options.bedfile)

    print(len(gc.chr_list))
    if len(gc.chr_list) == 1:
        gc = gc.chr_list[0]
    elif options.chromosome < 0 or options.chromosome > len(gc.chr_list):
        raise ValueError("invalid --chromosome value ; must be in [1-%s]" % len(gc.chr_list)+1)
    else:
        gc = gc.chr_list[options.chromosome-1]

    gc.running_median(n=options.ws)
    gc.coverage_scaling()
    gc.compute_zscore()
    gc.plot_coverage()
    from pylab import show
    show()


if __name__ == "__main__":
   import sys
   main()#sys.argv)

