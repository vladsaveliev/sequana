import argparse



class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_report"):
        usage = """Welcome to SEQUANA - Coverage standalone

            sequana_report --directory report

AUTHORS: Thomas Cokelaer, Dimitri Desvillechabrol
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("-d", "--directory", dest="directory", type=str,
            required=False, help="""filename of a BED file""")



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


if __name__ == "__main__":
   import sys
   main()#sys.argv)

