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

from sequana.reporting.report_main import BaseReport

from reports import HTMLTable


class SubMappingReport(BaseReport):
    """
    """
    def __init__(self, start, stop, chrom_index, chrom_name, thresholds,
            high_roi, low_roi, directory="report",
            output_filename="submapping.html", **kargs):
        super(SubMappingReport, self).__init__(
                jinja_filename="submapping/index.html",
                directory=directory, init_report=False,
                output_filename=output_filename, **kargs)
        self.jinja['title'] = "Mapping Report of {0} [{1},{2}]".format(
                chrom_name, start, stop)
        self.jinja['path'] = "../"
        self.chrom_index = chrom_index
        self.high_roi = high_roi
        self.low_roi = low_roi
        self.start = start
        self.stop = stop
        self.thresholds = thresholds

    def set_data(self, data):
        self.mapping = data

    def _get_region(self, df):
        return df[(df["end"] > self.start) & (df["start"] < self.stop)]

    def parse(self):
        self.mapping.write_csv(self.directory + os.sep +
                "mapping_{0}_{1}.{2}.csv".format(self.start, self.stop,
                self.chrom_index), start=self.start, stop=self.stop,
                header=False)
        self.jinja['input_df'] = "'mapping_{0}_{1}.{2}.csv'".format(self.start,
                self.stop, self.chrom_index)

        merge_low_cov = self._get_region(self.low_roi)
        self.jinja["low_cov_threshold"] = -self.thresholds.low
        self.jinja["low_cov_threshold_2"] = -self.thresholds.low2
        self.jinja["nb_low_region"] = len(merge_low_cov)
        html = HTMLTable(merge_low_cov)
        html.add_bgcolor("size")
        self.jinja['low_coverage'] = html.to_html(index=False)

        merge_high_cov = self._get_region(self.high_roi)
        self.jinja["high_cov_threshold"] = self.thresholds.high
        self.jinja["high_cov_threshold_2"] = self.thresholds.high2
        self.jinja["nb_high_region"] = len(merge_high_cov)
        html = HTMLTable(merge_high_cov)
        html.add_bgcolor("size")
        self.jinja['high_coverage'] = html.to_html(index=False)
