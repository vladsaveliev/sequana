# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>, 
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Main SEQUANA pipeline to create snakefile and executable scripts"""
import os
import shutil
import glob
import sys
import time
from optparse import OptionParser
import argparse

from easydev.console import red, purple, green, blue
from easydev import DevTools

import sequana
from sequana.snaketools import FastQFactory
from sequana.adapters import FindAdaptersFromIndex



class Tools(object):
    # Helper class to simplify following code
    dv = DevTools()
    def __init__(self, verbose=True):
        self.verbose = verbose
    def error(self, txt):
        print(red(txt))
        sys.exit(1)
    def purple(self, txt):
        if self.verbose: print(purple(txt))
    def red(self, txt):
        if self.verbose: print(red(txt))
    def green(self, txt):
        if self.verbose: print(green(txt))
    def blue(self, txt):
        if self.verbose: print(blue(txt))
    def onweb(self, link):
        from easydev import onweb
        onweb(link)
    def mkdir(self, name):
        self.dv.mkdir(name)
    def print(self, txt):
        if self.verbose: print(txt)


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

        group = self.add_argument_group('GENERAL')

        group.add_argument("--version", dest='version',
                action="store_true", help="print version")
        group.add_argument("--quiet", dest="verbose", default=True,
                          action="store_false")
        group.add_argument("--force-init", dest='force_init',
                action="store_true", 
                help="""Does not ask for permission to create the files (overwrite 
                    existing ones)""")
        group.add_argument("--issue", dest='issue',
                          action="store_true", help="Open github issue page")
        # Here we can use either --check-config <name> or  justg --check-config
        # So we must use nargs="?" but in such case, the default is stored in
        # const= so that if one does not use the --check-config, the default is
        # None but if one use it (without argument), the default is config.yaml
        group.add_argument("--check-config", dest='check_config',
                default=None, nargs="?", const="config.yaml" ,
                help="Check config file (e.g. that files exist or YAML is correct)")


        group = self.add_argument_group("PIPELINES")
        group.add_argument("--pipeline", dest='pipeline', type=str,
                required=False, 
                help="""Initialise a pipeline by fetching the corresponding 
                    snakefile and config file. Required files may also 
                    be copied in the project directory. Valid Pipelines' names 
                    can be obtained using --show-pipelines. See also --info 
                    <pipeline name> for more information about a specific 
                    pipeline. """)
        group.add_argument("--info", dest='info', type=str,
                required=False, 
                help="""Given a valid Sequana pipeline, this option opens the 
                    the corresponding README on GitHub""")
        group.add_argument("--show-pipelines", dest='show_pipelines',
                action="store_true", help="print names of the available pipelines")



        group = self.add_argument_group("CONFIGURATION")
        # options to fill the config file
        group.add_argument("--file1", dest="file1", type=str, 
            help=""" Fills the *samples:file1* field in the config file. To be used
                with --init option""")
        group.add_argument("--file2", dest="file2", type=str,
            help=""" Fills the *samples:file2* field in the config file. To be used
                with --init option""")
        group.add_argument("--glob", dest="glob", type=str,
                          help="a glob to find files. You can use wildcards")
        group.add_argument("--project", dest="project", type=str,
            help=""" Fills the *project* field in the config file. To be used
                with --init option""")
        group.add_argument("--kraken", dest="kraken", type=str,
            help=""" Fills the *kraken* field in the config file. To be used
                with --init option""")
        group.add_argument("--reference", dest="reference", type=str, 
            help=""" Fills the *reference* in bwa_ref section. To be used with --init option""")
        group.add_argument("--input-dir", dest="input_dir", type=str,
            help="""Search for a pair (or single) of reads in the directory, and
                fills automatically the project and file1/file2 fields in the config 
                file. To be used with --init option. If more than two files or not 
                files ending in fastq.gz are found, an error is raised.""")


        group = self.add_argument_group("ADAPTER and QUALITY related")
        group.add_argument("--no-adapters", dest="no_adapters", default=False,
            action="store_true")
        group.add_argument("--adapter-fwd", dest="adapter_fwd", type=str, 
            help="""A string representing the forward adapter. Can be a file 
                in FASTA format""")
        group.add_argument("--adapter-rev", dest="adapter_rev", type=str, 
            help="""A string representing the forward adapter. Can be a file 
                in FASTA format""")
        group.add_argument("--index-mapper", dest="index_mapper", type=str,
            help="""a CSV file with 3 columns named 'sample name', 'index1','index2' """)
        group.add_argument("--adapters", dest="adapters", type=str,
            help="""set to universal_nextera, universal_pcrfree""")
        group.add_argument("--quality-cutoff", dest="quality_cutoff", type=str,
                          default="30",
            help="""cutoff for single read or paired-end. If single read,
            provide a number between 0 and 40 e.g., 30 means remove all reads
            with quality below 30 (te details about the algorithm can be found
            in sequana documentation. For paired-end, provide a string as 30,30
            with two numbers separated by a comma (no space)""")
        group.add_argument("--quality", dest="quality_cutoff", type=str,
                          default="30",
                          help="""cutoff for single read or paired-end. If single read,
            provide a number between 0 and 40 e.g., 30 means remove all reads
            with quality below 30 (te details about the algorithm can be found
            in sequana documentation. For paired-end, provide a string as 30,30
            with two numbers separated by a comma (no space)""")

        # ====================================================== CLUSTER
        group = self.add_argument_group("Cluster", 
            """cluster related (LSF, slurm, ...)
            """)
        group.add_argument("--cluster", dest="cluster", type=str,
                          help="""a valid snakemake option dedicated to a cluster.
                          e.g on LSF cluster --cluster 'qsub -cwd -q<QUEUE> '"""
                          )
        group.add_argument("--jobs", dest="jobs", type=int, default=4,
                          help="""jobs to use on a cluster"""
                          )
        group.add_argument("--forceall", dest="forceall", action='store_true',
                          help="""add --forceall in the snakemake command"""
                          )


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

    sa = Tools(verbose=options.verbose)

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
    
    if options.info: 
        from sequana import Module
        module = Module(options.info)
        module.onweb()
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

    # In all other cases we must have either --pipeline, --run or --info (mutually
    # exclusive
    if options.pipeline and options.info:
        sa.error("ERROR: --pipeline and --info options are mutually exclusive")

    if options.pipeline:
        # check validity of the pipeline name
        if options.pipeline not in valid_pipelines:
            txt = "".join([" - %s\n" % this for this in valid_pipelines])
            sa.error("%s not a valid pipeline name. Use of one:\n" %
options.pipeline
                     + txt)

    # If user provides file1 and/or file2, check the files exist
    if options.file1 and os.path.exists(options.file1) is False:
        raise ValueError("%s does not exist" % options.file1)
    if options.file2 and os.path.exists(options.file2) is False:
        raise ValueError("%s does not exist" % options.file2)
    if options.kraken and os.path.exists(options.kraken) is False:
        raise ValueError("%s does not exist" % options.kraken)
    if options.adapter_rev and os.path.exists(options.adapter_rev) is False:
        sa.error('Invalid filename provided with --adapter-rev (must exists)')
    if options.adapter_fwd and os.path.exists(options.adapter_fwd) is False:
        sa.error('Invalid filename provided with --adapter-fwd (must exists)')


    if options.input_dir:
        options.glob = options.input_dir + os.sep + "*fastq.gz"


    # If a glob is used, we may have multiple project to handle
    if options.glob:
        ff = FastQFactory(options.glob)
        if options.index_mapper:
            adapter_finder = FindAdaptersFromIndex(options.index_mapper)
        elif options.adapter_fwd:
            pass
        elif options.adapters:
            pass
        elif options.no_adapters is True:
            pass
        else:
            sa.error("adapters need to be provided (or use --no-adapters")

        with open("multirun.sh", "w") as fout:
            import sequana
            fout.write("#!/usr/sh\n")
            fout.write("# generated with sequana version %s with this command\n" % sequana.version)
            fout.write("# %s\n" % " ".join(sys.argv))
            for tag in ff.tags:
                sa.print("Found %s project" % tag)
                options.project = tag
                options.file1 = ff.get_file1(tag)
                options.file2 = ff.get_file2(tag)
                if options.index_mapper:
                    fwd, rev = adapter_finder.save_adapters_to_csv(tag)
                    options.adapter_fwd = fwd
                    options.adapter_rev = rev
                sequana_init(options)
                fout.write("cd %s\n" % tag)
                fout.write("sh runme.sh &\n")
                fout.write("cd ..\n")
                fout.write("echo Starting %s\n" % tag)
                fout.write("sleep 0.5\n")
        if options.no_adapters is True and options.pipeline in ['quality', 'quality_taxon']:
            print("You did not provide information about adapters. You will have"
                "to edit the config.yaml file to fill that information")
        #Run sequana_init N times changing the files and project each time
        return

    # If --input-dir or --glob were not used, we now try to use the --file1 and --file2
    # options. We need to get the project name if not provided
    if options.file1 is None and options.file2 is None:
        if options.force_init is True:
            pass
        else:
            sa.error("You must provide one or two input files using"
                " --input-dir or --glob or a combo of --file1 and --file2. In "
                "the later case, you must provide at least --file1 (single-end)")
    elif options.file1 and options.file2 is None:
        ff = FastQFactory([options.file1])
        options.project = ff.tags[0]
    elif options.file1 and options.file2:
        ff = FastQFactory([options.file1, options.file2])
        options.project = ff.tags[0]

    if options.pipeline:
        sequana_init(options)
    elif options.info: 
        from sequana import Module
        module = Module(options.info)
        module.onweb()
        return

    if options.no_adapters is True and options.pipeline in ['quality', 'quality_taxon']:
        print("You did not provide information about adapters. You will have"
            "to edit the config.yaml file to fill that information")

def sequana_init(options):
    sa = Tools(verbose=options.verbose)
    import sequana.snaketools as sm
    from sequana import Module, SequanaConfig, sequana_data

    if options.project is None:
        options.project = input(red("No project name provided (use --project). "
            "Enter a project name:"))

    target_dir = options.project
    sa.verbose = True
    sa.blue("Creating project '" + target_dir +"'")
    sa.verbose = options.verbose

    try:
        module = Module(options.pipeline)
    except:
        sa.red("Invalid module name provided (%s). " % options.pipeline)
        sa.print("Valid pipeline names are:")
        for this in sm.pipeline_names:
            sa.print(" - %s" % this)
        sys.exit(1)

    # the module exists, so let us now copy the relevant files that is
    # the Snakefile, the config and readme files:
    if os.path.exists(target_dir):
        sa.print("Will override the following files if present: %s.rules" %
options.pipeline)
        if options.force_init is True:
            choice = "y"
        else:
            choice = input(red("Do you want to proceed ? [y]/n:"))
        if choice == "n":
            sys.exit(0)
    else:
        sa.print(purple("Creating %s directory" % options.project))
        sa.mkdir(options.project)

    sa.print("Copying snakefile")
    shutil.copy(module.snakefile, target_dir + os.sep + options.pipeline + ".rules")


    txt = "User command::\n\n"
    txt += "    %s \n\n" % " ".join(sys.argv)
    txt += "You can now run snakemake yourself or type::"
    txt += purple("""

    snakemake -s %s.rules --stats stats.txt -p -j 4

    """ % options.pipeline)
    txt += """
    # -j 4 means you will use 4 cores
    # -p prints the commands used
    # --stats stats.txt must be used since stats.txt is expected to be found.

    or just run the bash script::

        sh runme.sh


    EDIT THE config.yaml FILE TO SPECIFIED THE INPUT FILE LOCATION
    """

    sa.print("Creating README")
    with open(target_dir + os.sep + "README", "w") as fh:
        fh.write(txt.replace("\x1b[35m","").replace("\x1b[39;49;00m", ""))

    sa.print("Creating the config file")
    # Create (if needed) and update the config file
    config_filename = target_dir + os.sep + "config.yaml"
    #if os.path.exists(config_filename) is False or options.force_init:
    #    shutil.copy(module.config, target_dir + os.sep + "config.yaml")
    shutil.copy(module.config, target_dir + os.sep + "config.yaml")

    # Update the config file if possible
    with open(config_filename, "r") as fin:
        config_txt = fin.read()

    #print(options)
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
        else:
            params['kraken'] = None

        if options.reference:
            params['reference'] = os.path.abspath(options.reference)
        else:
            params['reference'] = None

        if options.file1 and options.adapter_fwd:
            params["adapter_fwd"] = "file:" + options.adapter_fwd
            os.rename(options.adapter_fwd, target_dir + os.sep + options.adapter_fwd)

        if options.file2  and options.adapter_rev:
            params["adapter_rev"] = "file:" + options.adapter_rev
            os.rename(options.adapter_rev, target_dir + os.sep + options.adapter_rev)

        if options.adapters == "universal":
            params["adapter_fwd"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC"
            params["adapter_rev"] = "TCTAGCCTTCTCGCAGCACATCCCTTTCTCACATCTAGAGCCACCAGCGGCATAGTAA"

        if options.quality_cutoff:
            # TODO handle single read
            if options.file2 is None:
                params['quality_cutoff'] = options.quality_cutoff
            else:
                params['quality_cutoff'] = "%s,%s" % (options.quality_cutoff, options.quality_cutoff)

        fout.write(config_txt % params)

    # figure out from the config file if any files are required
    try:
        # stop here if the file does not exists
        cfg = SequanaConfig(config_filename)
        if 'requirements' in cfg.config.keys():
            if "sequana" in cfg.config.requirements.keys():
                from sequana import sequana_data
                for filename in cfg.config.requirements.sequana:
                    sa.print('Copying %s from sequana' % filename)
                    shutil.copy(sequana_data(filename, "data"), target_dir)
            if "url" in cfg.config.requirements.keys():
                for link in cfg.config.requirements.url:
                    sa.print("This file %s will be needed" % filename)
    except Exception as err:
        print(err)
        print("No config file found")

    with open(target_dir + os.sep + "runme.sh", "w") as fout:
        import sequana
        cmd = "#!/usr/sh\n"
        cmd += "# generated with sequana version %s with this command\n"
        cmd += "# %s\n" % sys.argv
        cmd += "snakemake -s %(project)s.rules --stats report/stats.txt -p -j %(jobs)s"
        if options.forceall:
            cmd += " --forceall "

        if options.cluster:
            cmd += ' --cluster "%s"' % options.cluster

        cmd += " 1>run.out 2>run.err"
        fout.write(cmd % {'project':options.pipeline , 'jobs':options.jobs, 
			"version": sequana.version})

    sa.green("Initialisation of %s succeeded" % target_dir)
    sa.green("Please, go to the project directory ")
    sa.purple("\n   cd %s\n" % target_dir)
    sa.green("Check out the README and config.yaml files")
    sa.green("A basic script to run the analysis is named runme.sh ")
    sa.purple("\n    sh runme.sh\n")
    sa.green("In case of trouble, please post an issue on https://github.com/sequana/sequana/issue ")
    sa.green("or type sequana --issue and fill a post with the error and the config file (NO DATA PLEASE)")


    ## --cluster "qsub -cwd -qQUEUE -V -e -o "



if __name__ == "__main__":
   import sys
   main(sys.argv)

