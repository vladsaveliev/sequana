import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse



class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_summary"):
        usage = """Welcome to SEQUANA - Coverage standalone

            sequana_summary --bed file.bed --window-size 1001

AUTHORS: Thomas Cokelaer, Dimitri Desvillechabrol
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("-1", "--file1", dest="file1", type=str,
            required=True, help="""filename of a BED file""")
        self.add_argument("-2", "--file2", dest="file2", type=str,
            required=False, help="""filename of a BED file""")
        self.add_argument("-n", "--sample", default=500000, type=int)


def get_stats(filename, sample=500000):
    from sequana import FastQC
    ff = FastQC(filename, max_sample=sample)
    stats = ff.get_stats()
    return stats


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
    if options.file1 and not options.file2:
        stats = get_stats(options.file1, options.sample)
        print()
        for name in stats.columns:
            print("%s: %s" % (name, stats[name].values[0]))
    elif options.file1 and options.file2:
        from easydev import MultiProcessing
        mc = MultiProcessing(2)
        mc.add_job(get_stats, options.file1, options.sample)
        mc.add_job(get_stats, options.file2, options.sample)
        mc.run()
        df = mc.results[0]
        df = df.append(mc.results[1])

        print()
        print(df)


if __name__ == "__main__":
   import sys
   main()#sys.argv)

