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
import json
from .report_main import BaseReport

# externals
from easydev import precision
from reports import HTMLTable

import pandas as pd



class PhixReport(BaseReport):
    """Report dedicated to inform amount of phix in the FastQ


    """
    def __init__(self, output_filename="phix.html", directory="report",
            overwrite=False, **kargs):
        """

        :param jinja_template: name of a directory (either local) or
            from sequana/share/templates where JINJA files are available. A file
            named index.html is required but may be renamed (with
            **output_filename** parameter).
        :param output_filename: name of the final HTML file.
        :param directory: name of the output directory (defaults to report)

        Parameters accepted by :class:`reports.Report` are also accepted.

        """
        super(PhixReport, self).__init__(jinja_filename="phix_contaminant/index.html", 
                 directory=directory, output_filename=output_filename, **kargs)

        self.title = "Phix Report Summary"
        self.jinja['title'] = "Phix Report Summary"
        self.jinja['main_link'] = "index.html"
        self.input_filename = None # to be provided by the user

    def parse(self):

        if self.input_filename is None:
            return
        data = json.load(open(self.input_filename, "r"))
        for key, value in data.items():
            self.jinja[key] = value

        # Overwrite mode is a real name (rather than se or pe)
        if data['mode'] == "pe":
            self.jinja['mode'] = "Paired-end"
        elif data['mode'] == "se":
            self.jinja['mode'] = "Single-end"

        x = data['R1_mapped']
        y = data['R1_unmapped']

        # ad contamination inside jinja
        contamination = x / float(x+y) * 100
        self.jinja['contamination'] = precision(contamination, 3)

        # add HTML table 
        if "R2_mapped" in data.keys():
            df = pd.DataFrame({
                'R1': [data['R1_mapped'], data['R1_unmapped']],
                'R2': [data['R2_mapped'], data['R2_unmapped']]})
        else:
            df = pd.DataFrame({
                'R1': [data['R1_mapped'], data['R1_unmapped']]})
        df.index = ['mapped', 'unmapped']

        print(df)

        h = HTMLTable(df)
        html = h.to_html(index=True)

        html += "Unpaired: %s <hr>" % data['unpaired']
        html += "duplicated: %s <hr>" % data['duplicated']

        self.jinja['stats'] = html







