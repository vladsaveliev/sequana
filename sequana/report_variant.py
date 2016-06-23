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

import pandas as pd

from reports import HTMLTable
from sequana.report_main import BaseReport


class VariantReport(BaseReport):
    """
    """
    def __init__(self, csv_file, output_filename, directory="report", **kargs):
        super(VariantReport, self).__init__(jinja_filename="variant/index.html",
                directory=directory, output_filename=output_filename, 
                **kargs)
        self.jinja['title'] = "Variants Report"
        self.csv = csv_file

    def set_data(self, data):
        self.vcf_record = data

    def parse(self):
        df = self.vcf_record.df

        # Create readable CSV
        self.vcf_record.to_csv(self.csv)
        self.jinja["csv_link"] = self.csv.split("/")[-1]

        self.jinja["nb_variant"] = len(df)

        shutil.copy(self.vcf_record.filename, self.directory)
        self.jinja["vcf_link"] = self.vcf_record.filename.split("/")[-1] 

        html = HTMLTable(df)
        self.jinja['vcf_filter'] = html.to_html(index=False)
