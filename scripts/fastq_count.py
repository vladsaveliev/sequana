from sequana.fastq import FastQ

import sys
from optparse import OptionParser
import argparse


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="fastq_extract"):
        usage = """%s input N output \n""" % prog
        usage += """usage2: %s fastq_filename""" % prog
        usage += """Examples:

            fastq_extract input.fastq.gz
            fastq_extract input.fastq.gz 10000 output.fastq

        you can also use named arguments::

            fastq_extract --input input.fastq.gz --N 10000 --ouput output.fastq.gz

        """
        super(Options, self).__init__(usage=usage, prog=prog)
        self.add_argument("--input", dest='input_filename', type=str,
                            required=True, help="input fastq gzipped or not")
 
def main(args=None):
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="fastq_count")
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    elif len(args) == 2:
        class SimpleOpt():
            pass
        options = SimpleOpt()
        options.input_filename = args[1]
    else:
        options = user_options.parse_args(args[1:])

    f = FastQ(options.input_filename)
    # could be simplified calling count_reads only once
    print("Number of reads: %s" % f.count_reads())
    print("Number of lines %s " % f.count_lines())





if __name__ == "__main__":
   import sys
   main(sys.argv)
   
