import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse
from easydev.console import red, purple, green
import time

from easydev import DevTools


class Tools(object):
    dv = DevTools()
    def error(self, txt):
        print(red(txt))
        sys.exit(1)
    def purple(self, txt):
        print(purple(txt))
    def red(self, txt):
        print(red(txt))
    def green(self, txt):
        print(green(txt))
    def onweb(self, link):
        from easydev import onweb
        onweb(link)
    def mkdir(self, name):
        self.dv.mkdir(name)
sa = Tools()


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana"):
        usage = """Welcome to SEQUANA standalone

            sequana --init phix_removal
            sequana --init <sequana pipeline>  --file1 A.fastq.gz --project test
            sequana --init <sequana pipeline>  --input-dir ../
            sequana --show-pipelines
            sequana --version

        If you have lots of paired-end, use a --glob without --project ; This
        will create as many directories as required.

        Note that all files must have _R1_ and _R2_ tag in their names.
        
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
        self.add_argument("--issue", dest='issue',
                          action="store_true", help="Open github issue page")
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
        self.add_argument("--glob", dest="glob", type=str,
                          help="a glob to find files. You can use wildcards")
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

    # OPTIONS that gives info and exit

    if options.issue:
        sa.onweb('https://github.com/sequana/sequana//issues')
    if options.version:
        import sequana
        sa.purple("Sequana version %s" % sequana.version)
        return

    if options.show_pipelines:
        sa.purple("Valid pipeline names:")
        for this in sorted(valid_pipelines):
            print(" - " + this)
        return

    if options.check_config:
        from sequana import SequanaConfig
        cfg = SequanaConfig(options.check_config)
        # file1 must always be defined
        if os.path.exists(cfg.config.samples.file1) is False:
            print(red("%s does not exists. Please edit %s (samples:file1)" %
(cfg.config.samples.file1, options.check_config)))
            return

        sa.purple("The %s looks good" % options.check_config)
        return


    # In all other cases we must have either --init, --run or --info (mutually
    # exclusive
    if options.init and options.info:
        sa.error("ERROR: --init and --info options are mutually exclusive")


    if options.init:
        # check validity of the pipeline name
        if options.init not in valid_pipelines:
            txt = "".join([" - %s\n" % this for this in valid_pipelines])
            sa.error("%s not a valid pipeline name. Use of one:\n" % options.init
                     + txt)


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
        #if options.project is None:
        #    options.project = os.path.basename(os.path.abspath(options.input_dir))
        from sequana.snaketools import FastQFactory
        ff = FastQFactory(options.input_dir + os.sep + "*fastq.gz")
        options.project = ff.tags[0]
        files = glob.glob(options.input_dir + os.sep + "*.fastq.gz")
        assert len(files) <= 2
        if len(files) == 1:
            options.file1 = files[0]
            options.file2 = ""
            print("Found 1 fastq file %s" % options.file1)
        elif len(files) == 2:
            if "R1" in files[0]:
                options.file1 = files[0]
                options.file2 = files[1]
            elif "R2" in files[0]:
                options.file2 = files[0]
                options.file1 = files[1]
            print("Found 2 fastq files:")
            print("  - %s" % options.file1)
            print("  - %s" % options.file2)
        elif len(files) == 0:
            raise ValueError("No fastq.gz files found in %s " % options.input_dir)
        elif len(files) > 2:
            raise NotImplementedError("Ambiguous number of fastq greater than 2")

    if options.glob:
        from sequana.snaketools import FastQFactory
        ff = FastQFactory(options.glob)
        for tag in ff.tags:
            options.project = tag
            options.file1 = ff.get_file1(tag)
            options.file2 = ff.get_file2(tag)

            sequana_init(options)
        #Run sequana_init N times changing the files and project each time
        return

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


    if options.project is None:
        options.project = input(red("No project name provided (use --project). Enter a project name:"))

    target_dir = options.project


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
    if os.path.exists(target_dir):
        print("Will override the following files if present: %s.rules" % options.init)
        if options.force_init is True:
            choice = "y"
        else:
            choice = input(red("Do you want to proceed ? [y]/n:"))
        if choice == "n":
            sys.exit(0)
    else:
        print(purple("Creating %s directory" % options.project))
        sa.mkdir(options.project)

    time.sleep(0.5)
    print("Copying snakefile")
    shutil.copy(module.snakefile, target_dir + os.sep + options.init + ".rules")



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

        sh runme.sh


    EDIT THE config.yaml FILE TO SPECIFIED THE INPUT FILE LOCATION
    """

    print("Creating README")
    time.sleep(0.5)
    with open(target_dir + os.sep + "README", "w") as fh:
        fh.write(txt.replace("\x1b[35m","").replace("\x1b[39;49;00m", ""))


    print("Creating the config file")
    time.sleep(0.5)
    # Create (if needed) and update the config file
    config_filename = target_dir + os.sep + "config.yaml"
    if os.path.exists(config_filename) is False:
        shutil.copy(module.config, target_dir + os.sep + "config.yaml")

    # Update the config file if possible
    with open(config_filename, "r") as fin:
        config_txt = fin.read()

    with open(config_filename, "w") as fout:
        from collections import defaultdict
        params = defaultdict(str)
        params['project'] = options.project

        if options.file1:
            params['file1'] = os.path.abspath(options.file1)
        else:
            sa.green("You have not provided any input file. use --file1 or --input-dir")
            sa.green("You will need to edit the config.yaml file")
        if options.file2:
            params['file2'] = os.path.abspath(options.file2)
        else:
            sa.green("You have not provided any input file. use --file2 or --input-dir")
            sa.green("You will need to edit the config.yaml file")

        if options.kraken:
            params['kraken'] = os.path.abspath(options.kraken)
        if options.reference:
            params['reference'] = os.path.abspath(options.reference)

        fout.write(config_txt % params)

    # figure out from the config file if any files are required
    try:
        # stop here if the file does not exists
        cfg = SequanaConfig(config_filename)
        if 'requirements' in cfg.config.keys():
            if "sequana" in cfg.config.requirements.keys():
                from sequana import sequana_data
                for filename in cfg.config.requirements.sequana:
                    print('Copying %s from sequana' % filename)
                    shutil.copy(sequana_data(filename, "data"), target_dir)
            if "url" in cfg.config.requirements.keys():
                for link in cfg.config.requirements.url:
                    print("This file %s will be needed" % filename)
    except Exception as err:
        print(err)
        print("No config file found")

    with open(target_dir + os.sep + "runme.sh", "w") as fout:
        fout.write("#!/usr/sh\nsnakemake -s %s.rules --stats stats.txt -p -j 4" % options.init )

    sa.green("Initialisation of %s succeeded" % target_dir)
    sa.green("Please, go to the project directory ")
    sa.purple("\n   cd %s\n" % target_dir)
    sa.green("Check out the README and config.yaml files")
    sa.green("A basic script to run the analysis is named runme.sh ")
    sa.purple("\n    sh runme.sh\n")
    sa.green("In case of trouble, please post an issue on https://github.com/sequana/sequana/issue ")
    sa.green("or type sequana --issue and fill a post with the error and the config file (NO DATA PLEASE)")




if __name__ == "__main__":
   import sys
   main(sys.argv)

