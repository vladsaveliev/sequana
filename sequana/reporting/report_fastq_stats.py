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
import glob

from sequana.reporting.report_main import BaseReport

from sequana.lazy import reports
from sequana.lazy import pandas as pd


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
        self.filenames = []


    def parse(self):
        # Add the canvas js in the header
        from sequana.resources.canvas import bar

        if len(self.filenames) == 0:
            files = glob.glob("%s/*json" % self.input_directory)
        else:
            files = self.filenames

        acgt = [
            {"name":"A", "data": {}},
            {"name":"C", "data": {}},
            {"name":"G", "data": {}},
            {"name":"T", "data": {}}]

        dfsum = None
        for filename in files:
            try:
                # if the mapped file is empty, the file is empty
                thisdata = pd.read_json(filename)
            except:
                thisdata= pd.DataFrame(
                    {'A': {0: 0}, 'C': {0: 0}, 'G': {0: 0}, 'GC content': {0: 0},
                     'N': {0: 0}, 'T': {0: 0}, 'average read length': {0: 0},
                     'mean quality': {0: 0},
                     'n_reads': {0: 0},
                     'total bases': {0: 0}})

            if "R2_.unmapped" in filename:
                key = "R2.unmapped"
            elif "R1_.unmapped" in filename:
                key = "R1.unmapped"
            elif "R1_.mapped" in filename:
                key = "R1.mapped"
            elif "R2_.mapped" in filename:
                key = "R2.mapped"
            elif "_R1_" in filename:
                key = "R1"
            elif "_R2_" in filename:
                key = "R2"
            elif "_R1." in filename:
                key = "R1"
            elif "_R2." in filename:
                key = "R2"
            else:
                key = "R1"
            thisdata.index = [key]

            if dfsum is None:
                dfsum = thisdata.copy()
            else:
                dfsum = dfsum.append(thisdata)

            acgt[0]["data"][key] = thisdata['A'].values[0]
            acgt[1]["data"][key] = thisdata['C'].values[0]
            acgt[2]["data"][key] = thisdata['G'].values[0]
            acgt[3]["data"][key] = thisdata['T'].values[0]

        # First argument is just an Identifier
        data = bar.stacked_bar("Test", "ACGT distribution", acgt)

        data += """<script type="text/javascript" src="js/canvasjs.min.js"></script> """
        self.jinja["canvas"] = data



        dfsum["n_reads"] = [str(int(x)) for x in dfsum["n_reads"].values]
        dfsum["GC content"] = [str(int(x*100)/100.) for x in dfsum["GC content"].values]
        dfsum = dfsum.T
        dfsum = dfsum[sorted(dfsum.columns)]
        dfsum = dfsum.ix[['A', 'C', 'G', 'T', 'N', 'GC content', 'n_reads']]
        S = dfsum.ix[['A', 'C', 'G', 'T', 'N']].sum() 
        # FIXME change the json files themselves instead of multiplying by 100
        try:
            if S.values[0]>0:
                dfsum.ix[['A', 'C', 'G', 'T', 'N']] /= (S/100.) 
        except: 
            print("fixme in fastq_stats")

        html =""

        htmltable = reports.HTMLTable(dfsum).to_html(index=True)
        # Save the table into a temporary file
        with open(self.input_directory + "/temp.html", "w") as fh:
            fh.write(htmltable)

        html += htmltable
        html += """
    <div id="chartContainerTest" style="height: 300px; width: 100%;"></div>
        """

        # Copying images into the report/images directory
        #targets = self.copy_images_to_report("%s/images/*png" % self.input_directory)

        # Create table with those images
        #from sequana.htmltools import galleria
        #if len(targets):
        #    html += galleria(targets)
        self.jinja["content"] = html







