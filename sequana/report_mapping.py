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

"""Report dedicated to Mapping

"""
import os

import pandas as pd

from reports import HTMLTable
from sequana.report_main import BaseReport


class MappingReport(BaseReport):
    """Report dedicated to Mapping

    """
    def __init__(self, directory="report", project="",
            output_filename="report_mapping.html", **kargs):
        super(MappingReport, self).__init__(
                jinja_filename="mapping/index.html", directory=directory,
                output_filename=output_filename, **kargs)
        self.project = project
        self.jinja['title'] = "Mapping Report of {0}".format(project)

    def set_data(self, data):
        self.chr_list = data

    def parse(self):
        self.jinja['main_link'] = 'index.html'
        df = pd.DataFrame()
        formatter = '<a target="_blank" alt={0} href="{1}">{0}</a>'
        for chrom in self.chr_list:
            link = self.project + "_" + chrom.chrom_name + "_mapping.html"
            df = df.append({"chromosome": formatter.format(
                chrom.chrom_name, link), "size": "{0:,}".format(len(chrom))}, 
                ignore_index=True)
        html = HTMLTable(df)
        self.jinja['list_chromosome'] = html.to_html(index=False)
