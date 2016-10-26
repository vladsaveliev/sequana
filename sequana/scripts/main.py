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
from sequana.adapters import FindAdaptersFromIndex, AdapterReader
import sequana.snaketools as sm
from sequana import SequanaConfig, sequana_data


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



class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('FORMAT|'):
            return text[7:].splitlines() 
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana"):
        usage = """Welcome to SEQUANA standalone

        If you want to analyse several samples distributed in sub-directories,
        use:

            sequana --pipeline quality_control --pattern "./*/*fastq.gz" 
                --design SampleSheetUsed.csv  --adapters PCRFree

        If all samples are in the current directory:

            sequana --pipeline quality_control --input-directory . 
                --design SampleSheetUsed.csv  --adapters PCRFree

        If you want to analyse a specific pair of files:

            sequana --pipeline quality_control --file1 test_R1_.fastq.gz
                --file2 test_R2_.fastq.gz
                --design SampleSheetUsed.csv  --adapters PCRFree

        Note that files must end in fastq.gz, must contain _R1_ and _R2_ tag and
        must have a common sample name before _R1_ and _R2_ (like in the example
        above). You may have test after the _R1_ or _R2_ tag.

        Here above, --pipeline can refer to other pipelines such as
        variant_calling. Valid pipeline names can be retrieved using :

            sequana --show-pipelines

        If multiple samples are found, sub-directories are created for each
        sample. 

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
                description=description, formatter_class=SmartFormatter)

        group = self.add_argument_group('GENERAL')

        group.add_argument("-v", "--version", dest='version',
                action="store_true", help="print version")
        group.add_argument("-q", "--quiet", dest="verbose", default=True,
                          action="store_false")
        group.add_argument("-f", "--force", dest='force',
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

        group.add_argument("--get-config", dest="get_config", default=False,
                action="store_true",
                help=("Get config of a pipeline and copy it in the current "
                      "directory"))
        group.add_argument("--config-params", dest="config_params",
            type=str,
            help="""FORMAT|Overwrite any field in the config file by using
the following convention. A config file is in YAML format
and has a hierarchy of parametesr. For example:

    samples:
        file1: R1.fastq.gz
        file2: R2.fastq.gz
    bwa_mem_phix:
        mem:
            threads: 2

Here we have 3 sections with 1,2,3 levels respectively. On the
command line, each level is separated by a : sign and each
meter to be changed separated by a comma. So to change the
project name and threads inside the bwa_phix section use:

--config-params project:newname, bwa_mem_phix:mem:threads:4

Be aware that when using --config-params, all comments are removed.""")

        # ====================================================== CLUSTER
        group = self.add_argument_group("SNAKEMAKE AND CLUSTER RELATED",
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

        group = self.add_argument_group("INPUT FILES")
        group.add_argument("-1", "--file1", dest="file1", type=str,
            help=""" Fills the *samples:file1* field in the config file. To be used
                with --init option""")
        group.add_argument("-2", "--file2", dest="file2", type=str,
            help=""" Fills the *samples:file2* field in the config file. To be used
                with --init option""")
        group.add_argument("--pattern", dest="pattern", type=str,
            default="*.fastq.gz",
            help="a pattern to find files. You can use wildcards")
        group.add_argument("-i", "--input-directory", dest="input_directory", type=str,
            default=None,
            help="""Search for a pair (or single) of reads in the directory, and
                fills automatically the project and file1/file2 fields in the config
                file. To be used with --init option. If more than two files or not
                files ending in fastq.gz are found, an error is raised.""")
        group.add_argument("-o", "--output-directory", dest="target_dir", type=str,
            default="analysis",
            help="directory where to create files and store results")

        group.add_argument_group("SPECIFIC")
        group.add_argument("--reference", dest="reference", type=str,
            help=""" Fills the *reference* in bwa_ref section. To be used with --init option""")

        group = self.add_argument_group("QUALITY_CONTROL Pipeline")
        group.add_argument("--adapter-fwd", dest="adapter_fwd", type=str,
            help="""A string representing the forward adapter. Can be a file
                in FASTA format""")
        group.add_argument("--adapter-rev", dest="adapter_rev", type=str,
            help="""A string representing the forward adapter. Can be a file
                in FASTA format""")
        group.add_argument("--design", dest="design", type=str,
            help="""a CSV file with 3 columns named 'sample name', 'index1','index2' """)
        group.add_argument("--adapters", dest="adapters", type=str, default="",
            help="""When using --design, you must also provide the type of
                adapters using this parameter. Valid values are either Nextera or PCRFree
                Corresponding files can be found in github.com/sequana/sequana/resources/data/adapters
                """)
        group.add_argument("--kraken", dest="kraken", type=str,
            help=""" Fills the *kraken* field in the config file. To be used
                with --init option""")

def main(args=None):
    """Mostly checking the options provided by the user and then call
    :func:`sequana_init` function to create the pre-filled config file +
    snakemake + README +runme.sh in a dedicated project directory.

    """
    # these imports must be local
    from sequana.misc import textwrap
    from sequana.snaketools import Module

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        sa = Tools()
        sa.purple("Welcome to Sequana standalone application")
        sa.error("You must use --pipeline <valid pipeline name>\nuse "
                 "--show-pipelines or --help for more information")
        return
    else:
        options = user_options.parse_args(args[1:])

    sa = Tools(verbose=options.verbose)
    sa.purple("Welcome to Sequana standalone application")

    # We put the import here to make the --help faster
    from sequana.snaketools import pipeline_names as valid_pipelines

    # Those options are mutually exclusive
    flag = int("%s%s%s%s%s%s" % (
            int(bool(options.issue)),
            int(bool(options.version)),
            int(bool(options.info)),
            int(bool(options.show_pipelines)),
            int(bool(options.pipeline)),
            int(bool(options.get_config))
            ), 2)
    if flag not in [1,2,4,8,16,3]:
        sa.error("You must use one of --pipeline, --info, "
            "--show-pipelines, --issue, --version, --get-config")

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
            m = Module(this)
            sa.green(" - " + this)
            print(textwrap(m.overview, indent=8))
        return

    if options.info:
        module = Module(options.info)
        module.onweb()
        return

    if options.pipeline:
        # check validity of the pipeline name
        if options.pipeline not in valid_pipelines:
            txt = "".join([" - %s\n" % this for this in valid_pipelines])
            sa.error("%s not a valid pipeline name. Use of one:\n" % options.pipeline
                     + txt)

    # copy locally the request config file from a specific pipeline
    if flag == 3: #--get-config and --pipeline used
        module = Module(options.pipeline)
        copy_config_from_sequana(module)
        return

    # pipeline should be defined by now. Let us start the real work here
    Module("dag").check("warning")
    Module(options.pipeline).check("warning")

    # If user provides file1 and/or file2, check the files exist
    if options.file1 and os.path.exists(options.file1) is False:
        raise ValueError("%s does not exist" % options.file1)

    if options.file2 and os.path.exists(options.file2) is False:
        raise ValueError("%s does not exist" % options.file2)

    if options.kraken and os.path.exists(options.kraken) is False:
        raise ValueError("%s does not exist" % options.kraken)

    # check valid combo of arguments
    flag = int("%s%s%s%s%s" % (
            int(bool(options.pattern)),
            int(bool(options.input_directory)),
            int(bool(options.file1)),
            int(bool(options.file2)),
            int(bool(options.config)),
            ), 2)

    # config file has flag 1, others have flag 2,4,8
    # config file can be used with 2,4,8 so we also add 2+1, 4+1, 8+1
    if flag not in [1, 2, 4, 8, 3, 5, 9, 16,17,24,25]:
        sa.error(help_input + "\n\nUse --help for more information")


    # input_dir is nothing else than a glob to fastq.gz files
    if options.input_directory is None:
        options.input_directory = "./"
    options.glob = options.input_directory + os.sep + options.pattern

    def _get_adap(filename):
        return sequana_data(filename, "data/adapters")


    if options.verbose:
        print("Looking for sample files in %s" % options.glob)
    ff = FastQFactory(options.glob)
    if options.verbose:
        print("Found %s projects/samples " % len(ff.tags))

    if options.pipeline == 'quality_control':
        # check combo
        flag = int("%s%s%s%s" % (
            int(bool(options.design)),
            int(bool(options.adapters)),
            int(bool(options.adapter_fwd)),
            int(bool(options.adapter_rev))
            ), 2)

        if flag not in [12, 4, 2, 3]:
            sa.error("You must use a design experimental file using --design"
                     " and --adapters to indicate the type of adapters (PCRFree"
                     " or Nextera), or provide the adapters directly as a "
                     " string or a file using --adapter_fwd (AND --adapter_rev"
                     " for paired-end data). A third way is to set --adapters"
                     " to either Nextera, PCRFree or universal in which case "
                    " all adapters will be used (slower)")

        # flag 12 (design + adapters when wrong args provided)
        if options.design and options.adapters not in ["Nextera", "PCRFree"]:
            raise ValueError("When using --design, you must also "
                "provide the type of adapters using --adapters (set to "
                "Nextera or PCRFree)")
        # flag 12 (design + adapters with correct args)
        elif options.design and options.adapters in ["Nextera", "PCRFree"]:
            options.adapters_fwd = options.adapters
            options.adapters_rev = options.adapters
        else:
            if options.adapter_fwd is None:
                assert options.adapters in ["universal", "PCRFree", "Nextera"]
                # flag 4
                if options.adapters == "universal":
                    options.adapter_fwd = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC"
                    options.adapter_rev = "TCTAGCCTTCTCGCAGCACATCCCTTTCTCACATCTAGAGCCACCAGCGGCATAGTAA"
                # flag 4
                elif options.adapters == "PCRFree":
                    options.adapter_fwd = "file:" + _get_adap('adapters_PCR-free_PF1_220616_fwd.fa')
                    options.adapter_rev = "file:" + _get_adap('adapters_PCR-free_PF1_220616_rev.fa')
                # flag 4
                elif options.adapters == "Nextera":
                    options.adapter_fwd = "file:" + _get_adap("adapters_Nextera_PF1_220616_fwd.fa")
                    options.adapter_rev = "file:" + _get_adap("adapters_Nextera_PF1_220616_rev.fa")
            # flag 2/3
            else:
                if options.adapter_fwd:
                    # Could be a string or a file. If a file, check the format
                    if os.path.exists(options.adapter_fwd):
                        AdapterReader(options.adapter_fwd)
                        options.adapter_fwd = "file:%s" % options.adapter_fwd
                if options.adapter_rev:
                    # Could be a string or a file. If a file, check the format
                    if os.path.exists(options.adapter_rev):
                        AdapterReader(options.adapter_rev)
                        options.adapter_rev = "file:%s" % options.adapter_rev
        if options.design:
            # Just check the format
            adapter_finder = FindAdaptersFromIndex(options.design,
                                options.adapters)

    # If all options are valid, we can now create the tree structure
    sequana_init(options)


def copy_config_from_sequana(module, source="config.yaml", target="config.yaml"):
    # identify config name from the requested module
    user_config = module.path + os.sep + source
    if os.path.exists(user_config):
        shutil.copy(user_config, target)
        txt = "copied %s from sequana %s pipeline" % (source, module.name)
    else:
        txt = "%s does not exists locally or within Sequana "+\
              "library (%s pipeline)"
    print(txt)


def sequana_init(options):
    from sequana.misc import textwrap
    from sequana import Module
    sa = Tools(verbose=options.verbose)


    #sa.blue("Creating project directory (use --project to overwrite the " + \
    #        "inferred value): '" + options.project +"'", force=True)

    module = Module(options.pipeline)

    if os.path.exists(options.target_dir):
        txt = "Will override the following files if present: %s.rules " +\
              "config.yaml, runme.sh, ..."
        sa.blue(txt % options.pipeline)

        if options.force is True:
            choice = "y"
        else:
            choice = input(red("Do you want to proceed (to avoid this " +
                               " message, use --force)? [y]/n:"))

        if choice == "n":
            sys.exit(0)

    # Copying snakefile
    sa.print("Copying snakefile")
    sa.mkdir(options.target_dir)
    shutil.copy(module.snakefile, options.target_dir + os.sep + options.pipeline + ".rules")

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
    with open(options.target_dir + os.sep + "README", "w") as fh:
        fh.write(txt.replace("\x1b[35m","").replace("\x1b[39;49;00m", ""))

    # Creating Config file
    sa.print("Creating the config file")

    # Create (if needed) and update the config file
    config_filename = options.target_dir + os.sep + "config.yaml"

    if options.config:
        # full existing path
        if os.path.exists(options.config):
            shutil.copy(options.config, config_filename)
        else: # or a sequana config file in the module path ?
            copy_config_from_sequana(module, options.config, config_filename)
    else:
        shutil.copy(module.config, config_filename)

    # Update the config file if possible. first we read back the config file
    # requested
    with open(config_filename, "r") as fin:
        config_txt = fin.read()


    # and save it back filling it with relevant information provided by the user
    with open(config_filename, "w") as fout:
        from collections import defaultdict
        params = defaultdict(str)

        # TODO: copy design file in the working directory ?
        params['input_directory'] = os.path.abspath(options.input_directory)
        params['pattern'] = options.pattern
        if options.file1:
            params['file1'] = os.path.abspath(options.file1)
            params['input_directory'] = ""
        if options.file2:
            params['file2'] = os.path.abspath(options.file2)

        if options.pipeline == "quality_control":
            if options.design:
                shutil.copy(options.design, options.target_dir + os.sep )
                params['adapter_design'] = os.path.basename(options.design)
            else:
                params['adapter_design'] = ""

            if options.kraken:
                params['kraken'] = os.path.abspath(options.kraken)
            else:
                params['kraken'] = ""

            if options.reference:
                params['reference'] = os.path.abspath(options.reference)
            else:
                params['reference'] = None

            params["adapter_fwd"] = options.adapter_fwd
            params["adapter_rev"] = options.adapter_rev
            params["adapter_type"] = options.adapters

        fout.write(config_txt % params)

    # figure out from the config file if any files are required
    cfg = SequanaConfig(config_filename)
    if 'requirements' in cfg.config.keys():
        for requirement in cfg.config.requirements:
            if requirement.startswith('http') is False:
                from sequana import sequana_data
                sa.print('Copying %s from sequana' % requirement)
                shutil.copy(sequana_data(requirement, "data"), options.target_dir)
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
    with open(options.target_dir + os.sep + "runme.sh", "w") as fout:
        cmd = "#!/usr/sh\n"
        cmd += "# generated with sequana version %s with this command:\n" % sequana.version
        cmd += "# %s\n" % " ".join(sys.argv)
        cmd += "snakemake -s %(project)s.rules --stats stats.txt -p -j %(jobs)s --nolock"
        if options.forceall:
            cmd += " --forceall "

        if options.cluster:
            cmd += ' --cluster "%s"' % options.cluster

        if options.redirection:
            cmd += " 1>run.out 2>run.err"
        fout.write(cmd % {'project':options.pipeline , 'jobs':options.jobs,
            "version": sequana.version})

    sa.green("Initialisation of %s succeeded" % options.target_dir)
    sa.green("Please, go to the project directory ")
    sa.purple("\n   cd %s\n" % options.target_dir)
    sa.green("Check out the README and config.yaml files")
    sa.green("A basic script to run the analysis is named runme.sh ")
    sa.purple("\n    sh runme.sh\n")
    sa.green("In case of trouble, please post an issue on https://github.com/sequana/sequana/issue ")
    sa.green("or type sequana --issue and fill a post with the error and the config file (NO DATA PLEASE)")

    ## --cluster "qsub -cwd -qQUEUE -V -e -o "



if __name__ == "__main__":
   import sys
   main(sys.argv)

