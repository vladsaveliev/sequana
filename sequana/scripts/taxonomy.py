import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse





class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_taxonomy"):
        usage = """Welcome to SEQUANA - Taxonomy standalone

            sequana_taxonomy --file1 R1.fastq --file2 R2.fastq --database db

AUTHORS: Thomas Cokelaer
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("--file1", dest="file1", type=str,
            help="""R1 fastq file (zipped) """)
        self.add_argument("--file2", dest="file2", type=str,
            help="""R2 fastq file (zipped) """)
        self.add_argument("--database", dest="database", type=str,
            help="""Kraken DB """)
        self.add_argument("--output", dest="output", type=str,
            help="""name of the output HTML file""", default="krona.html")
        self.add_argument("--thread", dest="thread", type=int,
            help="""number of threads to use """, default=4)
        self.add_argument("--show-html", dest="html", 
            action="store_true", 
            help="""number of threads to use """, default=4)


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
    from sequana import KrakenPipeline

    fastq = []
    if options.file1:
        fastq.append(options.file1)
    if options.file2:
        fastq.append(options.file2)
    k = KrakenPipeline(fastq, options.database, threads=options.thread, 
            output=options.output)
    k.run()

    if options.html is True:
        k.show()

if __name__ == "__main__":
   import sys
   main(sys.argv)

