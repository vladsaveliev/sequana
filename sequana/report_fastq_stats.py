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
import easydev
import os
import glob
import json

from .report_main import BaseReport

# a utility from external reports package
from reports import HTMLTable

import pandas as pd


class FastQStatsReport(BaseReport):
    """


    """
    def __init__(self, input_directory, output_filename="fastq_stats.html",
                 directory="report", **kargs):
        """

        :param jinja_template: name of a directory (either local) or
            from sequana/share/templates where JINJA files are available. A file
            named index.html is required but may be renamed (with
            **output_filename** parameter).
        :param output_filename: name of the final HTML file.
        :param directory: name of the output directory (defaults to report)

        Parameters accepted by :class:`reports.Report` are also accepted.

        """
        super(FastQStatsReport, self).__init__(
                jinja_filename="fastq_stats/index.html",
                directory=directory,
                output_filename=output_filename, **kargs)
        self.input_directory = input_directory
        self.jinja['title'] = "FastQ Stats Report"

    def parse(self):
        # Add the canvas js in the header
        from sequana.resources.canvas import bar

        files = glob.glob("%s/*json" % self.input_directory)
        acgt = [
            {"name":"A", "data": {}},
            {"name":"C", "data": {}},
            {"name":"G", "data": {}},
            {"name":"T", "data": {}}]

        dfsum = None
        for filename in files:
            try:
                # if the mapped file is empty, the file is empty
                thisdata = json.load(open(filename))
            except:
                thisdata= {"A":0, "C":0, "G":0, "T":0, "n_reads":0, "GC content":0}

            if "R2.unmapped" in filename:
                key = "R2.unmapped"
            elif "R1.unmapped" in filename:
                key = "R1.unmapped"
            elif "R1.mapped" in filename:
                key = "R1.mapped"
            elif "R2.mapped" in filename:
                key = "R2.mapped"
            elif "_R1_" in filename:
                key = "R1"
            elif "_R2_" in filename:
                key = "R2"
            elif "_R1." in filename:
                key = "R1"
            elif "_R2." in filename:
                key = "R2"

            if dfsum is None:
                dfsum = pd.Series(thisdata).to_frame()
                dfsum.columns = [key]
            else:
                dfsum[key] = pd.Series(thisdata)

            acgt[0]["data"][key] = thisdata['A']
            acgt[1]["data"][key] = thisdata['C']
            acgt[2]["data"][key] = thisdata['G']
            acgt[3]["data"][key] = thisdata['T']

        # First argument is just an Identifier
        data = bar.stacked_bar("Test", "ACGT distribution", acgt)

        data += """<script type="text/javascript" src="js/canvasjs.min.js"></script> """
        self.jinja["canvas"] = data

        dfsum = dfsum[sorted(dfsum.columns)]
        dfsum = dfsum.ix[['A', 'C', 'G', 'T', 'GC content', 'n_reads']]
        S = dfsum.ix[['A', 'C', 'G', 'T']].sum()
        dfsum.ix[['A', 'C', 'G', 'T']] /= S

        html = HTMLTable(dfsum).to_html(index=True)

        html += """
    <div id="chartContainerTest" style="height: 300px; width: 100%;"></div>
        """


        # Assuming in fastq directory, we figure out the HTML files and
        # create a table accordingly.
        sources = glob.glob("%s/*png" % self.input_directory)
        import shutil
        targets = []
        for source in sources:
            filename = source.rsplit("/", 1)[1]
            target = "report/images/" + filename
            shutil.copy(source, target)
            targets.append(target.replace("report/", ""))

        from sequana.htmltools import galleria
        html += galleria(targets)

        self.jinja["content"] = html







