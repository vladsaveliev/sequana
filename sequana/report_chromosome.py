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

"""Report dedicated to Mapping for one chromosome.

"""
import os

import pandas as pd

from reports import HTMLTable
from sequana.report_main import BaseReport
from sequana.report_submapping import SubMappingReport


class ChromosomeMappingReport(BaseReport):
    """Report dedicated to Mapping for one chromosome.

    """
    def __init__(self, chrom, low_threshold=-3, high_threshold=3, 
            directory="report", project="", **kargs):
        """.. rubric:: constructor

        """
        output_filename = project + "_" + chrom + "_mapping.html"
        super(ChromosomeMappingReport, self).__init__(
                jinja_filename="chromosome/index.html",
                directory=directory, 
                output_filename=output_filename, **kargs)
        self.project = project
        self.low_t = low_threshold
        self.high_t = high_threshold

    def set_data(self, data):
        self.mapping = data

    def _generate_submapping(self, low_df, high_df):
        i=0
        while True:
            i += 1
            if len(self.mapping) / i < 500000:
                step = int(len(self.mapping) / i) + 1
                break
        df = pd.DataFrame()
        formatter = '<a target="_blank" alt={0} href="{1}">{0}</a>'
        for i in range(0, len(self.mapping), step):
            name = "{0}_{1}_mapping_{2}".format(self.project, 
                    self.mapping.chrom_name, i)
            stop = i + step
            if stop > len(self.mapping):
                stop = len(self.mapping)
            name = "{0}_{1}".format(name, stop)
            link = name + ".html"
            r = SubMappingReport(start=i, stop=stop, low_df=low_df,
                    high_df=high_df, output_filename=name + ".html",
                    directory=self.directory, low_threshold=self.low_t,
                    high_threshold=self.high_t)
            r.jinja["main_link"] = "index.html"
            r.set_data(self.mapping)
            r.create_report()
            df = df.append({"name": formatter.format(name, link)}, 
                ignore_index=True)
        return df

    def parse(self):
        self.jinja['title'] = "Mapping Report of {0}".format(
                self.mapping.chrom_name)
        self.jinja['main_link'] = 'index.html'

        # Coverage plot
        self.jinja["cov_plot"] = self.project + "_" + self.mapping.chrom_name +\
                "_coverage.png"
        self.mapping.plot_coverage(filename=self.directory + os.sep + 
                self.jinja["cov_plot"], low_threshold=self.low_t, 
                high_threshold=self.high_t)

        # Barplot of zscore
        bp_paragraph = ("Distribution of the z-score (normalised coverage); "
                        "You should see a Gaussian distribution centered around 0"
            "The estimated parameters are mu={0:.2f} and sigma={1:.2f}.")
        self.jinja["bp_paragraph"] = bp_paragraph.format(
                self.mapping.best_gaussian["mu"], 
                self.mapping.best_gaussian["sigma"])
        self.jinja["bp_plot"] =  self.project + \
                "_" + self.mapping.chrom_name + "_zscore_hist.png"
        self.mapping.plot_hist(filename= self.directory + os.sep + 
                self.jinja["bp_plot"])

        # Low threshold case
        low_cov_paragraph = ("Regions with a z-score lower than {0:.2f} and at "
            "least one base with a z-score lower than {1:.2f} are detected as "
            "low coverage region. Thus, there are {2} low coverage regions")
        low_cov_df = self.mapping.get_low_coverage(self.low_t / 2)
        merge_low_cov = low_cov_df.merge_region(self.low_t)
        self.jinja["lc_paragraph"] = low_cov_paragraph.format(self.low_t / 2,
                self.low_t, len(merge_low_cov))
        html = HTMLTable(merge_low_cov)
        html.add_bgcolor("size")
        self.jinja['low_coverage'] = html.to_html(index=False)

        # High threshold case
        high_cov_paragraph = ("Regions with a z-score higher than {0:.2f} and at "
            "least one base with a z-score higher than {1:.2f} are detected as "
            "high coverage region. Thus, there are {2} high coverage regions")
        high_cov_df = self.mapping.get_high_coverage(self.high_t / 2)
        merge_high_cov = high_cov_df.merge_region(self.high_t)
        self.jinja['hc_paragraph'] = high_cov_paragraph.format(self.high_t / 2, 
                self.high_t, len(merge_high_cov))
        html = HTMLTable(merge_high_cov)
        html.add_bgcolor("size")
        self.jinja['high_coverage'] = html.to_html(index=False)

        # Sub mapping with javascript
        df = self._generate_submapping(merge_low_cov, merge_high_cov)
        html = HTMLTable(df)
        self.jinja['list_submapping'] = html.to_html(index=False) 
