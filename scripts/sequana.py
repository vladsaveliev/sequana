import os
import shutil
import sys
from optparse import OptionParser
import argparse


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana"):
        usage = """%s init fix_removal\n""" % prog
        usage += """Examples:

            sequana init <sequana pipeline>
            sequana run +snakemake options
            sequana report

        """
        super(Options, self).__init__(usage=usage, prog=prog)


        self.add_argument("--init", dest='snakefile', type=str,
                required=False,
                help="Get the snakefile and config file and possible other files")
        self.add_argument("--info", dest='info', type=str,
                required=False, help="Open README on the web")
        self.add_argument("--version", dest='version',
                action="store_true", help="print version")


def main(args=None):
    from sequana.snaketools import modules, Module, SequanaConfig
    import shutil
    from easydev import DevTools
    devtools = DevTools()

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
       options = user_options.parse_args(args[1:])

    if options.version:
        import sequana
        print(sequana.version)
        sys.exit(0)

    if options.info:
        module = Module(options.info)
        module.onweb()
        sys.exit(0)

    try:
        module = Module(options.snakefile)
    except:
        print("Invalid module name provided (%s). " % options.snakefile)
        print("Here are the valid modules from sequana: \n - %s" %
"\n - ".join(sorted(modules.keys())))
        sys.exit(1)

    # the module exists, let us now copy the relevant files that is
    # the Snakefile, the config and readme:
    print("Will override the following files: config.yaml, README.rst,"
          " Snakefile")
    from easydev.console import red
    choice = input(red("Do you want to proceed ? y/n"))
    if choice == "y":
        pass
    else:
        sys.exit(0)

    shutil.copy(module.snakefile, "Snakefile")
    shutil.copy(module.config, "config.yaml")
    shutil.copy(module.readme, "README.rst")

    # we also crate the directory to hold the raw data
    devtools.mkdir("fastq_raw")

    print("""You can now run snakemake yourself or type::

    snakemake -s Snakefile --stats stats.txt -p -j 4

    # -j 4 means you will use 4 cores
    # -p prints the commands used
    # --stats stats.txt must be used since stats.txt is expected to be found.

    or just run the bash script::

        sh sequana.sh


    EDIT THE config.yaml FILE TO SPECIFIED THE INPUT FILE LOCATION
    """)


    # figure out from the config file if any files are required
    try:
        # stop here if the file does not exists
        configfile = os.path.split(module.config)[1]
        config = SequanaConfig(configfile)
        if 'requirements' in config.parameters.keys():
            from sequana import sequana_data
            for filename in config.parameters.requirements:
                try:
                    print('Getting %s ' % filename)
                    shutil.copy(sequana_data(filename, "data"), ".")
                except:
                    print("This file %s will be needed" % filename)

    except Exception as err:
        print(err)
        print("No config file found")











if __name__ == "__main__":
   import sys
   main(sys.argv)

