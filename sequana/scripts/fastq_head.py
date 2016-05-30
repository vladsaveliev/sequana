from sequana.fastq import FastQ

from optparse import OptionParser
import argparse


class Options(argparse.ArgumentParser):
    def  __init__(self,  prog="fastq_extract"):
        usage = """%s input N output \n""" % prog
        usage += """usage2: %s --input input --N N --output output""" % prog
        usage += """Examples:

            fastq_extract input.fastq.gz 10000 output.fastq.gz
            fastq_extract input.fastq.gz 10000 output.fastq
            fastq_extract input.fastq 10000 output.fastq.gz
            fastq_extract input.fastq 10000 output.fastq

        you can also use named arguments::
           
            fastq_extract --input input.fastq.gz --N 10000 --ouput output.fastq.gz

        """
        super(Options, self).__init__(usage=usage,  prog=prog)
        self.add_argument("--nlines", dest='N', type=int, required=True,
                          help="Number of lines to extract.")
        self.add_argument("--input", dest='input_filename', type=str,
                            required=True, help="input fastq gzipped or not")
        self.add_argument("--output", dest='output_filename', type=str, 
                            required=True, 
                            help="output file with .gz extension or not")
 
def main(args=None):
    import sys
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="fastq_extract")
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    elif len(args) == 4:
        class SimpleOpt():
            pass
        options = SimpleOpt()
        options.input_filename = args[1]
        options.N = int(args[2])
        options.output_filename = args[3]
    else:
        options = user_options.parse_args(args[1:])

    f = FastQ(options.input_filename)
    f.extract_head(N=options.N,
                   output_filename=options.output_filename)





if __name__ == "__main__":
   import sys
   main(sys.argv)
   
