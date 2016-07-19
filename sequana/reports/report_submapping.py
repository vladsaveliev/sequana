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

from reports import HTMLTable
from sequana.report_main import BaseReport


class SubMappingReport(BaseReport):
    """
    """
    def __init__(self, start, stop, low_threshold, high_threshold, 
            low_df, high_df, directory="report", 
            output_filename="submapping.html", **kargs):
        super(SubMappingReport, self).__init__(
                jinja_filename="submapping/index.html",
                directory=directory,
                output_filename=output_filename, **kargs)
        self.jinja['title'] = "Mapping Report [{0},{1}]".format(start, stop)
        self.low_df = low_df
        self.high_df = high_df
        self.start = start
        self.stop = stop
        self.low_t = low_threshold
        self.high_t = high_threshold

    def set_data(self, data):
        self.mapping = data

    def _get_region(self, df):
        return df[(df["stop"] > self.start) & (df["start"] < self.stop)]

    def parse(self):
        self.mapping.write_csv(self.directory + os.sep +
                "mapping_{0}_{1}.csv".format(self.start, self.stop), 
                start=self.start, stop=self.stop, header=False)
        self.jinja['input_df'] = "'mapping_{0}_{1}.csv'".format(self.start, 
                self.stop)

        merge_low_cov = self._get_region(self.low_df)
        self.jinja["low_cov_threshold"] = self.low_t
        self.jinja["low_cov_threshold_2"] = "{0:.2f}".format(
                float(self.low_t) / 2)
        self.jinja["nb_low_region"] = len(merge_low_cov)
        html = HTMLTable(merge_low_cov)
        html.add_bgcolor("size")
        self.jinja['low_coverage'] = html.to_html(index=False)
        
        merge_high_cov = self._get_region(self.high_df)
        self.jinja["high_cov_threshold"] = self.high_t
        self.jinja["high_cov_threshold_2"] = "{0:.2f}".format(
                float(self.high_t) / 2)
        self.jinja["nb_high_region"] = len(merge_high_cov)
        html = HTMLTable(merge_high_cov)
        html.add_bgcolor("size")
        self.jinja['high_coverage'] = html.to_html(index=False)
