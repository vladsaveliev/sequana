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
import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_taxonomy"):
        usage = """Welcome to SEQUANA - Taxonomy standalone

            sequana_taxonomy --file1 R1.fastq --file2 R2.fastq --database db

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
            help="""Kraken DB """)
        self.add_argument("--output", dest="output", type=str,
            help="""name of the output HTML file""", default="kraken.html")
        self.add_argument("--thread", dest="thread", type=int,
            help="""number of threads to use """, default=4)
        self.add_argument("--show-html", dest="html", 
            action="store_true", 
            help="""number of threads to use """)


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
    from sequana import KrakenPipeline
    from easydev import DevTools
    devtools = DevTools()

    fastq = []
    if options.file1:
        fastq.append(options.file1)
    if options.file2:
        fastq.append(options.file2)


    devtools.mkdir("kraken")
    k = KrakenPipeline(fastq, options.database, threads=options.thread, 
           output="kraken" + os.sep + options.output)

    output_png = "kraken.png"
    k.run(output_png="kraken/%s" % output_png)

    if options.html is True:
        k.show()

    if 1==1:
        # Here we create a simple temporary config file to be read by the Summary
        # report 
        from easydev import TempFile
        config_txt = "samples:\n"
        config_txt += '    file1: "%s"\n' % options.file1
        if options.file2:
            config_txt += '    file2: "%s"\n'% options.file2
        config_txt += "project: Taxonomic Content\n"
        config_txt += "kraken:\n"
        config_txt += "    do: yes\n"
        config_txt += "    database: %s\n" % options.database

        tf = TempFile()
        fh = open(tf.name, "w")
        fh.write(config_txt)
        fh.close()
        print(tf.name)

        from sequana import SequanaSummary
        ss = SequanaSummary("kraken", "summary.html", tf.name, 
            include_all=False, workflow=False)
        ss.include_input_links()
        ss.jinja['kraken_pie'] = output_png
        ss.create_report()


if __name__ == "__main__":
   import sys
   main(sys.argv)

