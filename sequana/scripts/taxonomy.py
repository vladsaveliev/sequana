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
"""Standalone dedicated to taxonomic content (kraken)"""
import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse

from easydev import DevTools
from sequana import logger

class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_taxonomy"):
        usage = """Welcome to SEQUANA - Taxonomy standalone

        This standalone takes as input one or two (paired end) FastQ files.
        Using Kraken, Krona and some codecs from Sequana, it tries to identify
        the taxon/organism that match each reads found in the input files.

        Kraken requires an input database. We provide 3 databases but you can
        also use Sequana to help you building a customised one.

        We provide a DB called toydb. It contains only a few
        measles viruses. Its size is only 32Mb and should be used for testing
        and examples only.

        Another database is the so-called minikraken provided by Kraken's
        authors. It is about 4Gb and contains viruses and bacteria only.

        A third database is built with Sequana and is about 8Gb. It is
        stored on Synapse website and you will need an account (synapse.org).
        It contains about 22,000 species: viruses, bacteria, homo sapiens, fungi,

        Each DB can be downloaded using:

            sequana_taxonomy --download toydb

        Then, you need to use this kind of command:

            sequana_taxonomy --file1 R1.fastq --file2 R2.fastq
                --database /home/user/.config/sequana/kraken_toydb
                --show-html --thread 4

AUTHORS: Thomas Cokelaer
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:
        """
        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("--file1", dest="file1", type=str,
            help="""R1 fastq file (zipped) """)
        self.add_argument("--file2", dest="file2", type=str,
            help="""R2 fastq file (zipped) """)
        self.add_argument("--database", dest="database", type=str,
            #choices=["sequana_db1", "toydb", "minikraken"],
            help="""Path to a valid Kraken database. If you do not have any, use
                --download option""")
        self.add_argument("--output-directory", dest="directory", type=str,
            help="""name of the output directory""", default="taxonomy")
        self.add_argument("--thread", dest="thread", type=int,
            help="""number of threads to use """, default=4)
        self.add_argument("--show-html", dest="html",
            action="store_true",
            help="""Results are stored in report/ directory and results are
                not shown by default""")
        self.add_argument("--download", dest="download", type=str,
            default=None, choices=["sequana_db1", "toydb", "minikraken"],
            help="""download an official sequana DB. The sequana_db1 is stored
                in a dedicated Synapse page (www.synapse.org). minikraken
                is donwload from the kraken's author page, and toydb from
                sequana github.""")
        self.add_argument('--verbose', dest="verbose", default=False,
            action="store_true")


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
    from sequana import KrakenPipeline, SequanaConfig
    devtools = DevTools()

    if options.download:
        from sequana import KrakenDownload
        kd = KrakenDownload()
        kd.download(options.download)
        sys.exit()

    fastq = []
    if options.file1:
        devtools.check_exists(options.file1)
        fastq.append(options.file1)
    if options.file2:
        devtools.check_exists(options.file2)
        fastq.append(options.file2)


    from sequana import sequana_config_path as scfg
    if options.database == "toydb":
        options.database = "kraken_toydb"
    elif options.database == "minikraken":
        options.database = "minikraken_20141208"
    else:
        if options.database is None:
            logger.critical("You must provide a database")
            sys.exit(1)

    if os.path.exists(scfg + os.sep + options.database): # in Sequana path
        options.database = scfg + os.sep + options.database
    elif os.path.exists(options.database): # local database
        pass
    else:
        msg = "Invalid database name (%s). Neither found locally "
        msg += "or in the sequana path %s; Use the --download option"
        raise ValueError(msg % (options.database, scfg))

    output_directory = options.directory + os.sep + "kraken"
    devtools.mkdirs(output_directory)

    # if DB exists locally, use it otherwise add the sequana path
    k = KrakenPipeline(fastq, options.database, threads=options.thread,
        output_directory=output_directory, verbose=options.verbose)

    k.run()

    cfg = SequanaConfig()
    cfg.config.input_directory = None
    cfg.config.input_samples = {"file1": options.file1, "file2": options.file2}
    cfg.config.input_pattern = None
    cfg.config.input_extension = None
    cfg.config.kraken = {"database": options.database}

    from sequana.reporting.report_summary import SequanaSummary
    from sequana.snaketools import DummyManager

    filenames = [options.file1]
    if options.file2:
        filenames.append(options.file1)
    manager = DummyManager(filenames)

    sample = "custom"
    ss = SequanaSummary(sample, options.directory, output_filename="summary.html",
        configfile=cfg,
        include_all=False, workflow=False, manager=manager)

    ss.include_kraken()
    ss.create_report()

    if options.verbose:
        logger.info("Open ./%s/summary.html" % options.directory)
        logger.info("or ./%s/kraken/kraken.html" % options.directory)

    if options.html is True:
        ss.onweb()

if __name__ == "__main__":
   import sys
   main(sys.argv)

