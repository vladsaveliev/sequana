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
import sequana.snaketools as sm
from sequana import Module, SequanaConfig, sequana_data


help_input = """Incorrect combo of parameters.

You must provide one of those combinations:

  - one or two input files using --file1 (and --file2) if file2 provided, file1
    must be provided
  - the directory where to find the files with --input-dir
  - a pattern using --glob surrounded with quotes e.g. "*.fastq.gz"

An alternative is to use a local config file with the field file1 and file2
already filled. If you use --config with --file1 (and --file2), then the
pre-filled fields will be overwritten.

The input-dir is the same as --glob except that input-dir takes all gz whereas
--glob takes a more specific pattern e.g. *AB*fastq.gz

Note that --glob does not allow the --project to be used  (ignored)

This is also true for the --project that replaces project and other parameters
that are used to fill the config file:
"""


class Tools(object):
    # Helper class to simplify following code
    dv = DevTools()
    def __init__(self, verbose=True):
        self.verbose = verbose
    def error(self, txt):
        print(red(txt))
        sys.exit(1)
    def purple(self, txt, force=False):
        if self.verbose or force is True: 
            print(purple(txt))
    def red(self, txt, force=False):
        if self.verbose or force is True: 
            print(red(txt))
    def green(self, txt, force=False):
        if self.verbose or force is True: 
            print(green(txt))
    def blue(self, txt, force=False):
        if self.verbose or force is True: 
            print(blue(txt))
    def onweb(self, link):
        from easydev import onweb
        onweb(link)
    def mkdir(self, name):
        self.dv.mkdir(name)
    def print(self, txt, force=False):
        if self.verbose or force: print(txt)


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana"):
        usage = """Welcome to SEQUANA standalone

            sequana --pipeline quality
            sequana --pipeline quality  --file1 A.fastq.gz --project myname
            sequana --pipeline variant_calling  --input-dir ../
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
        group.add_argument("--force", dest='force',
                action="store_true",
                help="""Does not ask for permission to create the files (overwrite
                    existing ones)""")
        group.add_argument("--redirection", dest="redirection",
            default=False, action="store_true")

        group.add_argument("--issue", dest='issue',
                          action="store_true", help="Open github issue page")


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
        group.add_argument('--config', dest="config", type=str,
            default=None,
            help="""The default config file can be overwritten with a local file
(if found) or from a sequana file. For instance, if the quality_taxon module
directory contains a config file named config_test.yaml, then you can use it
with this option by just giving the name. Given the pipeline name, this utility
will fetch the config file automatically from sequana library.""")

        group.add_argument("--file1", dest="file1", type=str,
            help=""" Fills the *samples:file1* field in the config file. To be used
                with --init option""")
        group.add_argument("--file2", dest="file2", type=str,
            help=""" Fills the *samples:file2* field in the config file. To be used
                with --init option""")
        group.add_argument("--glob", dest="glob", type=str,
                          help="a glob to find files. You can use wildcards")
        group.add_argument("--input-dir", dest="input_dir", type=str,
            help="""Search for a pair (or single) of reads in the directory, and
                fills automatically the project and file1/file2 fields in the config
                file. To be used with --init option. If more than two files or not
                files ending in fastq.gz are found, an error is raised.""")

        group.add_argument("--project", dest="project", type=str,
            help=""" Fills the *project* field in the config file. To be used
                with --pipeline option""")
        group.add_argument("--kraken", dest="kraken", type=str,
            help=""" Fills the *kraken* field in the config file. To be used
                with --init option""")
        group.add_argument("--reference", dest="reference", type=str,
            help=""" Fills the *reference* in bwa_ref section. To be used with --init option""")

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
            help="""When using --index-mapper, you must also provide the type of 
                adapters using this parameter. Valid values are either Nextera or PCRFree
                Corresponding files can be found in github.com/sequana/sequana/resources/data/adapters
                """)

        group.add_argument("--config-params", dest="config_params", 
            type=str, 
            help="""Overwrite any field in the config file by using
                    the following convention. A config file is in YAML format
                    and has a hierarchy of parametesr. For example:

                    project: tutorial
                    samples:
                        file1: R1.fastq.gz
                        file2: R2.fastq.gz
                    bwa_phix:
                        mem:    
                            threads: 2

                Here we have 3 sections with 1,2,3 levels respectively. On the
                command line, each level is separated by a : sign and each
                meter to be changed separated by a comma. So to change the
                project name and threads inside the bwa_phix section use:

                --config-params project:newname, bwa_phix:mem:threads:4
                
                Be aware that when using --config-params, all comments are
removed.""")

        # ====================================================== CLUSTER
        group = self.add_argument_group("Snakemake and cluster related",
            """When launching the pipeline, one need to use the snakemake
command, which is written in the runme.sh file. One can tune that command to
define the type of cluster command to use
            """)
        group.add_argument("--snakemake-cluster", dest="cluster", type=str,
                          help="""a valid snakemake option dedicated to a cluster.
                          e.g on LSF cluster --cluster 'qsub -cwd -q<QUEUE> '"""
                          )
        group.add_argument("--snakemake-jobs", dest="jobs", type=int, default=4,
                          help="""jobs to use on a cluster"""
                          )
        group.add_argument("--snakemake-forceall", dest="forceall", action='store_true',
                          help="""add --forceall in the snakemake command"""
                          )


def main(args=None):
    """Mostly checking the options provided by the user and then call
    :func:`sequana_init` function to create the pre-filled config file +
    snakemake + README +runme.sh in a dedicated project directory.

    """
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        sa = Tools()
        sa.purple("Welcome to Sequana standalone application")
        sa.error("You must use --pipeline <valid pipeline name>\nuse --show-pipelines or --help for more information")
        return
    else:
        options = user_options.parse_args(args[1:])

    sa = Tools(verbose=options.verbose)
    sa.purple("Welcome to Sequana standalone application")

    # We put the import here to make the --help faster
    from sequana.snaketools import pipeline_names as valid_pipelines

    # Those options are mutually exclusive
    flag = int("%s%s%s%s%s" % (
            int(bool(options.issue)),
            int(bool(options.version)),
            int(bool(options.info)),
            int(bool(options.show_pipelines)),
            int(bool(options.pipeline)),
            ), 2)
    if flag not in [1,2,4,8,16]:
        sa.error("You must use one of --pipeline, --info, "
            "--show-pipelines, --issue, --version ")

    # OPTIONS that gives info and exit
    if options.issue:
        sa.onweb('https://github.com/sequana/sequana/issues')

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

    # In all other cases we must have either --pipeline, --run or --info (mutually
    # exclusive
    if options.pipeline and options.info:
        sa.error("ERROR: --pipeline and --info options are mutually exclusive")

    if options.pipeline:
        # check validity of the pipeline name
        if options.pipeline not in valid_pipelines:
            txt = "".join([" - %s\n" % this for this in valid_pipelines])
            sa.error("%s not a valid pipeline name. Use of one:\n" % options.pipeline
                     + txt)

    # If user provides file1 and/or file2, check the files exist
    if options.file1 and os.path.exists(options.file1) is False:
        raise ValueError("%s does not exist" % options.file1)

    if options.file2 and os.path.exists(options.file2) is False:
        raise ValueError("%s does not exist" % options.file2)

    if options.kraken and os.path.exists(options.kraken) is False:
        raise ValueError("%s does not exist" % options.kraken)

    if options.adapter_rev and os.path.exists(options.adapter_rev) is False:
        with open("adapter_rev.fa", "w") as fout:
            fout.write(">user\n%s" % options.adapter_rev)
        options.adapter_rev = "adapter_rev.fa"
    if options.adapter_fwd and os.path.exists(options.adapter_fwd) is False:
        with open("adapter_fwd.fa", "w") as fout:
            fout.write(">user\n%s" % options.adapter_fwd)
        options.adapter_fwd = "adapter_fwd.fa"

    # check valid combo of --glob / --fileX --input-dir
    # the 3 options are mutually exclusive
    flag = int("%s%s%s%s" % (
            int(bool(options.input_dir)),
            int(bool(options.glob)),
            int(bool(options.file1)),
            int(bool(options.config)),
            ), 2)

    # config file has flag 1, others have flag 2,4,8
    # config file can be used with 2,4,8 so we also add 2+1, 4+1, 8+1
    if flag not in [1, 2, 4, 8, 3, 5, 9]:
        sa.error(help_input + "\n\nUse --help for more information")

    if options.project is None:
        sa.blue("Note that --project was not provided. Will infer project " +\
                "name from the prefix of the filenames (before first " +\
                 "underscore)")

    # input_dir is nothing else than a glob to fastq.gz files
    if options.input_dir:
        options.glob = options.input_dir + os.sep + "*fastq.gz"

    # If a glob is used, we may have multiple project to handle
    # This is done using the FastQFactory class
    if options.glob:
        ff = FastQFactory(options.glob)
        if options.index_mapper:
            if options.adapters is None or options.adapters not in ["Nextera", "PCRFree"]:
                raise ValueError("When using --index-mapper, you must also " 
                    "provide the type of adapters using --adapters (set to "
                    "Nextera or PCRFree)")
            adapter_finder = FindAdaptersFromIndex(options.index_mapper, 
                                options.adapters)
        elif options.adapter_fwd:
            pass
        elif options.adapters:
            pass
        elif options.no_adapters is True:
            pass

        with open("multirun.sh", "w") as fout:
            import sequana
            fout.write("#!/usr/sh\n")
            fout.write("# generated with sequana version %s with this command\n" % sequana.version)
            fout.write("# %s\n" % " ".join(sys.argv))
            for tag in ff.tags:
                sa.print("Found %s project" % tag)
                #if options.project is None:
                options.project = tag
                options.file1 = ff.get_file1(tag)
                options.file2 = ff.get_file2(tag)
                if options.index_mapper:
                    fwd, rev = adapter_finder.save_adapters_to_fasta(tag)
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
    if options.file1 is None and options.file2 is None and options.config is None:
        sa.error(help_input)
    elif options.file1 and options.file2 is None:
        ff = FastQFactory([options.file1])
        if options.project is None:
            options.project = ff.tags[0]
    elif options.file1 and options.file2:
        ff = FastQFactory([options.file1, options.file2])
        if options.project is None:
            options.project = ff.tags[0]

    if options.pipeline:
        sequana_init(options)

    if options.no_adapters is True and options.pipeline in ['quality', 'quality_taxon']:
        sa.red("You did not provide information about adapters. You will have"
            "to edit the config.yaml file to fill that information")


def sequana_init(options):
    sa = Tools(verbose=options.verbose)

    if options.project is None:
        options.project = input(red("No project name provided and could not" + \
            " infer it (use --project). Enter a project name:"))

    target_dir = options.project
    sa.blue("Creating project directory (use --project to overwrite the " + \
            "inferred value): '" + options.project +"'", force=True)

    module = Module(options.pipeline)

    # the module exists, so let us now copy the relevant files that is
    # the Snakefile, the config and readme files:
    if os.path.exists(options.project):
        txt = "Will override the following files if present: %s.rules " +\
              "config.yaml, runme.sh, ..."
        sa.blue(txt % options.pipeline)

        if options.force is True:
            choice = "y"
        else:
            choice = input(red("Do you want to proceed ? [y]/n:"))

        if choice == "n":
            sys.exit(0)
    else:
        sa.purple("Creating %s directory" % options.project)
        sa.mkdir(options.project)

    # Copying snakefile
    sa.print("Copying snakefile")
    shutil.copy(module.snakefile, target_dir + os.sep + options.pipeline + ".rules")

    # Creating README to print on the screen and in a file
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

    EDIT THE config.yaml if needed

    Once finished with success, the report/ directory contains a summary.html
    and relevant files (depends on the pipeline).
    """

    sa.print("Creating README")
    with open(target_dir + os.sep + "README", "w") as fh:
        fh.write(txt.replace("\x1b[35m","").replace("\x1b[39;49;00m", ""))

    # Creating Config file
    sa.print("Creating the config file")

    # Create (if needed) and update the config file
    config_filename = target_dir + os.sep + "config.yaml"

    if options.config:
        if os.path.exists(options.config):
            shutil.copy(options.config, config_filename)
        else:
            # identify config name from the requested module
            user_config = module.path + os.sep + options.config
            if os.path.exists(user_config):
                shutil.copy(user_config, config_filename)
            else:
                txt = "%s does not exists locally or within Sequana "+\
                      "library (%s pipeline)"
                raise FileExistsError(txt % (user_config, options.pipeline))
    else:
        shutil.copy(module.config, config_filename)

    # Update the config file if possible. first we read back the config file
    # requested
    with open(config_filename, "r") as fin:
        config_txt = fin.read()

    # and save it back filling it with relevant information such as 
    # - project
    # - file1 
    # - file2 
    # - kraken db
    # - reference for the bwa_ref
    # - adapter
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
            shutil.move(options.adapter_fwd, target_dir + os.sep + options.adapter_fwd)

        if options.file2  and options.adapter_rev:
            params["adapter_rev"] = "file:" + options.adapter_rev
            shutil.move(options.adapter_rev, target_dir + os.sep + options.adapter_rev)

        if options.adapters == "universal":
            params["adapter_fwd"] = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC"
            params["adapter_rev"] = "TCTAGCCTTCTCGCAGCACATCCCTTTCTCACATCTAGAGCCACCAGCGGCATAGTAA"

        fout.write(config_txt % params)

    # figure out from the config file if any files are required
    cfg = SequanaConfig(config_filename)
    if 'requirements' in cfg.config.keys():
        for requirement in cfg.config.requirements:
            if requirement.startswith('http') is False:
                from sequana import sequana_data
                sa.print('Copying %s from sequana' % requirement)
                shutil.copy(sequana_data(requirement, "data"), target_dir)
            elif requirement.startswith("http"):
                from sequana.misc import wget
                sa.print("This file %s will be needed" % requirement)
                wget(requirement)

    # FIXME If invalid, no error raised 
    if options.config_params:
        cfg = SequanaConfig(config_filename)
        params = [this.strip() for this in options.config_params.split(",")]
        for param in params:
            if param.count(":") not in [1,2, 3]:
                txt = "incorrect format following --config-params"
                txt += "Expected at least one : sign or at most 2 of them"
                txt += "Config file section such as :\n"
                txt += "project: tutorial\n"
                txt += "should be encoded project:tutorial"
                raise ValueError(txt)
            if param.count(":") == 1:
                k,v = param.split(':')
                cfg.config[k] = v
            elif param.count(":") == 2:
                k1,k2,v = param.split(":")
                cfg.config[k1][k2] = v
            elif param.count(":") == 3:
                k1,k2,k3,v = param.split(":")
                cfg.config[k1][k2][k3] = v
        cfg.save(config_filename)

    # Creating a unique runme.sh file
    with open(target_dir + os.sep + "runme.sh", "w") as fout:
        cmd = "#!/usr/sh\n"
        cmd += "# generated with sequana version %s with this command:\n" % sequana.version
        cmd += "# %s\n" % " ".join(sys.argv)
        cmd += "snakemake -s %(project)s.rules --stats report/stats.txt -p -j %(jobs)s"
        if options.forceall:
            cmd += " --forceall "

        if options.cluster:
            cmd += ' --cluster "%s"' % options.cluster

        if options.redirection:
            cmd += " 1>run.out 2>run.err"
        fout.write(cmd % {'project':options.pipeline , 'jobs':options.jobs,
            "version": sequana.version})


    with open(target_dir + os.sep + "cleanme.py", "w") as fout:
        fout.write("""
import glob
import os
import shutil
from easydev import shellcmd
import time

directories = glob.glob("*")

for this in directories:
    if os.path.isdir(this) and this not in ['logs', 'data', 'report']:
        print('Deleting %s' % this)
        time.sleep(0.2)
        shellcmd("rm -rf %s" % this)
shellcmd("rm -f *rules README runme.sh *fa config.yaml snakejob.*")
shellcmd("rm -f cleanme.py")
shellcmd("rm -f runme.sh.*")
shellcmd("rm -rf .snakemake")
""")

    sa.green("Initialisation of %s succeeded" % target_dir)
    sa.green("Please, go to the project directory ")
    sa.purple("\n   cd %s\n" % target_dir)
    sa.green("Check out the README and config.yaml files")
    sa.green("A basic script to run the analysis is named runme.sh ")
    sa.purple("\n    sh runme.sh\n")
    sa.green("In case of trouble, please post an issue on https://github.com/sequana/sequana/issue ")
    sa.green("or type sequana --issue and fill a post with the error and the config file (NO DATA PLEASE)")


    ## --cluster "qsub -cwd -qQUEUE -V -e -o "


def check_config(config):
    from sequana import SequanaConfig
    cfg = SequanaConfig(config)
    # file1 must always be defined
    if os.path.exists(cfg.config.samples.file1) is False:
        print(red("%s does not exists. Please edit %s (samples:file1)" %
(cfg.config.samples.file1, check_config)))
        return
    sa.purple("The %s looks good" % check_config)



if __name__ == "__main__":
   import sys
   main(sys.argv)

