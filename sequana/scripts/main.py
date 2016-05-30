import os
import shutil
import sys
from optparse import OptionParser
import argparse





class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana"):
        usage = """%s init fix_removal\n""" % prog
        usage += """Examples:

            sequana --init <sequana pipeline>
            sequana --run <sequana pipeline>.rules 
            sequana --show-pipelines
            sequana --version

        """
        super(Options, self).__init__(usage=usage, prog=prog)

        self.add_argument("--init", dest='init', type=str,
                required=False, 
                help="Get the snakefile and config file and possible other files")
        self.add_argument("--run", dest='run', type=str,
                required=False, 
                help="Get the snakefile and config file and possible other files")
        self.add_argument("--info", dest='info', type=str,
                required=False, help="Open README on the web")
        self.add_argument("--version", dest='version',
                action="store_true", help="print version")
        self.add_argument("--show-pipelines", dest='show_pipelines',
                action="store_true", help="print available pipelines")
        self.add_argument("--force-init", dest='force_init',
                action="store_true", help="print available pipelines")


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
    from sequana.snaketools import pipeline_names as valid_pipelines
    from easydev.console import red, purple

    # OPTIONS that gives info and exit

    if options.version:
        import sequana
        print(purple("Sequana version %s" % sequana.version))
        return


    if options.show_pipelines:
        print(purple("Valid pipeline names:"))
        for this in sorted(valid_pipelines):
            print(" - " + this)
        return


    # In all other cases we must have either --init, --run or --info (mutually
    # exclusive
    if options.init and options.info:
        print(red("ERROR: --init and --info options are mutually exclusive"))
        sys.exit(1)
    if options.init and options.run:
        print(red("ERROR: --init and --run options are mutually exclusive"))
        sys.exit(1)
    if options.run and options.info:
        print(red("ERROR: --run and --info options are mutually exclusive"))
        sys.exit(1)


    if options.init:
        sequana_init(options)
    elif options.run:
        raise NotImplementedError
    elif options.info: 
        from sequana import Module
        module = Module(options.info)
        module.onweb()
        return


def sequana_init(options):

    import sequana.snaketools as sm
    from sequana import Module, SequanaConfig
    #modules, Module
    from easydev.console import red, purple

    try:
        module = Module(options.init)
    except:
        print(red("Invalid module name provided (%s). " % options.init))
        print("Valid pipeline names are:")
        for this in sm.pipeline_names:
            print(" - %s" % this)
        sys.exit(1)

    # the module exists, let us now copy the relevant files that is
    # the Snakefile, the config and readme:
    print("Will override the following files: config.yaml,"
          " %s.rules if they  exist" % options.init)

    if options.force_init is True:
        choice = "y"
    else:
        choice = input(red("Do you want to proceed ? y/n"))
        if choice == "y":
            pass
        else:
            sys.exit(0)

    shutil.copy(module.snakefile, options.init + ".rules")
    shutil.copy(module.config, "config.yaml")


    txt = "You can now run snakemake yourself or type::"
    txt += purple("""

    snakemake -s %s.rules --stats stats.txt -p -j 4

    """ % options.init)
    txt += """
    # -j 4 means you will use 4 cores
    # -p prints the commands used
    # --stats stats.txt must be used since stats.txt is expected to be found.

    or just run the bash script::

        sh sequana.sh


    EDIT THE config.yaml FILE TO SPECIFIED THE INPUT FILE LOCATION
    """
    print(txt)

    with open("README", "w") as fh:
        fh.write(txt.replace("\x1b[35m","").replace("\x1b[39;49;00m", ""))


    # figure out from the config file if any files are required
    try:
        # stop here if the file does not exists
        configfile = os.path.split(module.config)[1]
        cfg = SequanaConfig(configfile)
        if 'requirements' in cfg.config.keys():
            if "sequana" in cfg.config.requirements.keys():
                from sequana import sequana_data
                for filename in cfg.config.requirements.sequana:
                    print('Copying %s from sequana' % filename)
                    shutil.copy(sequana_data(filename, "data"), ".")
            if "url" in cfg.config.requirements.keys():
                for link in cfg.config.requirements.url:
                    print("This file %s will be needed" % filename)

    except Exception as err:
        print(err)
        print("No config file found")





if __name__ == "__main__":
   import sys
   main(sys.argv)

