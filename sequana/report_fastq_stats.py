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

        for filename in files:
            thisdata = json.load(open(filename))
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
            print(filename)
            acgt[0]["data"][key] = thisdata['A']
            acgt[1]["data"][key] = thisdata['C']
            acgt[2]["data"][key] = thisdata['G']
            acgt[3]["data"][key] = thisdata['T']

            # { "n_reads": 1627, "GC content": 0.45281800009750867}

        # First argument is just an Identifier
        data = bar.stacked_bar("Test", "ACGT distribution", acgt)

        data += """<script type="text/javascript" src="js/canvasjs.min.js"></script> """
        self.jinja["canvas"] = data


        html = """
    <div id="chartContainerTest" style="height: 300px; width: 100%;"></div>
        """


        # Assuming in fastq directory, we figure out the HTML files and
        # create a table accordingly.
        links = glob.glob("%s/*png" % self.input_directory)
        links = [x.split("/",1)[1] for x in links]

        from sequana.htmltools import galleria
        html += galleria(links)

        self.jinja["content"] = html







