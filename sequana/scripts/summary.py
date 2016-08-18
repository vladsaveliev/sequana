import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse



class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_summary"):
        usage = """Welcome to SEQUANA - Summary standalone

            sequana_summary --file file.fastq 
            sequana_summary --glob "file*.fastq"

AUTHORS: Thomas Cokelaer, Dimitri Desvillechabrol
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:

        prints basic stats about a set of input files
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("-f", "--file", dest="file", type=str,
            required=False, help="""filename of a FastQ file""")
        self.add_argument("-g", "--glob", dest="glob", type=str,
            required=False, help="""a glob/pattern of files. Must use quotes
                e.g. "*.fastq.gz" """)
        self.add_argument("-n", "--sample", default=500000, type=int)


def get_fastq_stats(filename, sample=500000):
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
    if options.file:
        options.glob = options.file

    from easydev import MultiProcessing
    mc = MultiProcessing(4)
    filenames = glob.glob(options.glob)
    for filename in filenames:
        mc.add_job(get_fastq_stats, filename, options.sample)
    mc.run()
    for i, this in enumerate(filenames):
        if i == 0:
            df = mc.results[0]
        else:
            df = df.append(mc.results[i])
    df.index = filenames

    print()
    print(df)


if __name__ == "__main__":
   import sys
   main()#sys.argv)

