import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse





class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana"):
        usage = """Welcome to SEQUANA standalone

            sequana --init phix_removal
            sequana --init <sequana pipeline>  --file1 A.fastq.gz --project test
            sequana --init <sequana pipeline>  --input-dir ../
            sequana --show-pipelines
            sequana --version
        
AUTHORS: Thomas Cokelaer and Dimitri Devilleschabrol
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION: the SEQUANA standalone provides tools to 
        fetch NGS pipelines developped within the SEQUANA porject. The main 
        option is the --init option that fetches a pipeline and the main config 
        file (config.yaml). The pipelines are SNAKEFILE and can be executed 
        with a specific command provided in the README (typically snakemake -s 
        pipeline.rules). For more information, you may use the --info followed
        by a pipeline name. You can also have a look at the online documentation
        on http://sequana.readthedocs.io or the source code on GitHub at  
        http://github.com/sequana/sequana
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        self.add_argument("--init", dest='init', type=str,
                required=False, 
                help="""Get the snakefile and config file corresponding to the
                    pipeline. Possibly other required files may be downloaded
                    Pipelines' names can be obtained using --show-pipelines
                    option. See also --info <pipeline name> for more information
                    about a specific pipeline. This option copies locally the
                    config.yam, README.rst and the requested pipelines 
                    <pipeline name>.rules""")
        self.add_argument("--info", dest='info', type=str,
                required=False, 
                help="""Given a known Sequana pipeline, this option opens the 
                    the corresponding README on GitHub""")
        self.add_argument("--show-pipelines", dest='show_pipelines',
                action="store_true", help="print names of the available pipelines")
        self.add_argument("--force-init", dest='force_init',
                action="store_true", 
                help="""Does not ask for permission to create the files (overwrite 
                    existing ones)""")
        self.add_argument("--version", dest='version',
                action="store_true", help="print version")

        # Here we can use either --check-config <name> or  justg --check-config
        # So we must use nargs="?" but in such case, the default is stored in
        # const= so that if one does not use the --check-config, the default is
        # None but if one use it (without argument), the default is config.yaml
        self.add_argument("--check-config", dest='check_config',
                default=None, nargs="?", const="config.yaml" ,
                help="Check config file (e.g. that files exist or YAML is correct)")


        # options to fill the config file
        self.add_argument("--file1", dest="file1", type=str, 
            help=""" Fills the *samples:file1* field in the config file. To be used
                with --init option""")
        self.add_argument("--file2", dest="file2", type=str,
            help=""" Fills the *samples:file2* field in the config file. To be used
                with --init option""")
        self.add_argument("--project", dest="project", type=str,
            help=""" Fills the *project* field in the config file. To be used
                with --init option""")
        self.add_argument("--kraken", dest="kraken", type=str,
            help=""" Fills the *kraken* field in the config file. To be used
                with --init option""")
        self.add_argument("--reference", dest="reference", type=str, 
            help=""" Fills the *reference* in bwa_ref section. To be used with --init option""")
        self.add_argument("--input-dir", dest="input_dir", type=str,
            help="""Search for a pair (or single) of reads in the directory, and
                fills automatically the project and file1/file2 fields in the config 
                file. To be used with --init option. If more than two files or not 
                files ending in fastq.gz are found, an error is raised.""")


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

    print(options.check_config)
    if options.check_config:
        from sequana import SequanaConfig
        cfg = SequanaConfig(options.check_config)
        # file1 must always be defined
        if os.path.exists(cfg.config.samples.file1) is False:
            print(red("%s does not exists. Please edit %s (samples:file1)" %
(cfg.config.samples.file1, options.check_config)))
            return

        print(purple("The %s looks good" % options.check_config))
        return


    # In all other cases we must have either --init, --run or --info (mutually
    # exclusive
    if options.init and options.info:
        print(red("ERROR: --init and --info options are mutually exclusive"))
        sys.exit(1)

    # If user provides file1 and/or file2, check the files exist
    if options.file1 and os.path.exists(options.file1) is False:
        raise ValueError("%s does not exist" % options.file1)
    if options.file2 and os.path.exists(options.file2) is False:
        raise ValueError("%s does not exist" % options.file2)
    if options.kraken and os.path.exists(options.kraken) is False:
        raise ValueError("%s does not exist" % options.kraken)
    if options.input_dir and os.path.isdir(options.input_dir) is False:
        raise ValueError("%s does not exist" % options.input_dir)
    if options.input_dir:
        # The directory should be used to fill automatically the files
        # and project name. The project name will be the directory name
        options.project = os.path.basename(os.path.abspath(options.input_dir))
        files = glob.glob(options.input_dir + os.sep + "*.fastq.gz")
        assert len(files) <= 2
        if len(files) == 1:
            options.file1 = files[0]
            options.file2 = ""
        elif len(files) == 2:
            if "R1" in files[0]:
                options.file1 = files[0]
                options.file2 = files[1]
            elif "R2" in files[0]:
                options.file2 = files[0]
                options.file1 = files[1]
        elif len(files) == 0 or len(files) >2:
            raise ValueError("No fastq.gz files found in %s " % options.input_dir)


    if options.init:
        sequana_init(options)
    elif options.info: 
        from sequana import Module
        module = Module(options.info)
        module.onweb()
        return


def sequana_init(options):

    import sequana.snaketools as sm
    from sequana import Module, SequanaConfig, sequana_data
    from easydev.console import red, purple

    try:
        module = Module(options.init)
    except:
        print(red("Invalid module name provided (%s). " % options.init))
        print("Valid pipeline names are:")
        for this in sm.pipeline_names:
            print(" - %s" % this)
        sys.exit(1)

    # the module exists, so let us now copy the relevant files that is
    # the Snakefile, the config and readme files:
    print("Will override the following files if present: config.yaml,"
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


    txt = "User command::\n\n"
    txt += "    %s \n\n" % " ".join(sys.argv)
    txt += "You can now run snakemake yourself or type::"
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


    fin = open("config.yaml", "r")
    config_txt = fin.read()
    fin.close()
    with open("config.yaml", "w") as fout:
        fout.write(config_txt % {
            'file1': options.file1,
            'file2': options.file2,
            'project': options.project,
            'kraken': options.kraken,
            'reference': options.reference
            }
        )

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

