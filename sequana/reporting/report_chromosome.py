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
from sequana.reporting.report_main import BaseReport
from sequana.reporting.report_submapping import SubMappingReport


class ChromosomeMappingReport(BaseReport):
    """Report dedicated to Mapping for one chromosome.

    """
    def __init__(self, mapping, features=None,
            directory="report", sample="", verbose=True, **kargs):
        """.. rubric:: constructor

        """
        self.mapping = mapping
        output_filename = "{0}_mapping.chrom{1}.html".format(sample, 
                mapping.chrom_index)
        super(ChromosomeMappingReport, self).__init__(
                jinja_filename="chromosome/index.html",
                directory=directory,
                output_filename=output_filename, **kargs)
        self.chrom_index = mapping.chrom_index
        self.sample = sample
        self.features = features
        self._max_window = 500000
        self.jinja['chrom_name'] = self.mapping.chrom_name
        self.verbose = verbose

    def _get_step(self):
        i=0
        while True:
            i += 1
            if len(self.mapping) / i < self._max_window:
                step = int(len(self.mapping) / i) + 1
                break

        return self._max_window
        return step

    def _generate_submapping(self, high_roi, low_roi):
        step = self._get_step()

        df = pd.DataFrame()
        formatter = '<a target="_blank" alt={0} href="{1}">{0}</a>'
        for i in range(0, len(self.mapping), step):
            if self.verbose:
                print("Generating sub mapping {}".format(i))
            name = "{0}_mapping_{1}".format(self.sample, i)
            stop = i + step
            if stop > len(self.mapping):
                stop = len(self.mapping)
            name = "{0}_{1}".format(name, stop)
            link = "submapping/{0}.chrom{1}.html".format(name, self.chrom_index)
            r = SubMappingReport(start=i, stop=stop, high_roi=high_roi, 
                    low_roi=low_roi, chrom_index=self.chrom_index,
                    chrom_name=self.mapping.chrom_name,
                    output_filename=name + ".chrom%i.html" % self.chrom_index,
                    directory=self.directory + "/submapping", 
                    thresholds=self.mapping.thresholds)
            r.jinja["main_link"] = self.filename
            r.set_data(self.mapping)
            r.create_report()
            df = df.append({"name": formatter.format(name, link)}, 
                ignore_index=True)
        return df

    def parse(self):
        self.jinja['title'] = "Mapping Report of {0}".format(
                self.mapping.chrom_name)
        self.jinja['nav_off'] = True

        # Stats of chromosome
        if self.verbose:
            print("Creating stats")
        df = pd.Series(self.mapping.get_stats()).to_frame()
        df.index.name = "Metric"
        df = df.reset_index()
        self.jinja["nc_stats"] = HTMLTable(df).to_html(index=False, header=True)

        # Coverage plot
        if self.verbose:
            print("Creating coverage plot")
        self.jinja["cov_plot"] = "images/{0}_coverage.chrom{1}.png".format(
                self.sample, self.chrom_index)
        self.mapping.plot_coverage(filename=self.directory + os.sep + 
                self.jinja["cov_plot"])

        # Barplot of normalized coverage with predicted gaussians
        if self.verbose:
            print("Creating zscore plots")
        nc_paragraph = ("Distribution of the normalized coverage with "
                "predicted Gaussian. The red line should be followed the "
                "trend of the barplot.")
        self.jinja["nc_paragraph"] = nc_paragraph
        self.jinja["nc_plot"] = "images/{0}_norm_cov_his.chrom{1}.png".format(
                self.sample, self.chrom_index)
        self.mapping.plot_hist_normalized_coverage(filename=self.directory +
                os.sep + self.jinja["nc_plot"])

        # Barplot of zscore
        bp_paragraph = ("Distribution of the z-score (normalised coverage); "
            "You should see a Gaussian distribution centered around 0. "
            "The estimated parameters are mu={0:.2f} and sigma={1:.2f}.")
        self.jinja["bp_paragraph"] = bp_paragraph.format(
                self.mapping.best_gaussian["mu"], 
                self.mapping.best_gaussian["sigma"])
        self.jinja["bp_plot"] = "images/{0}_zscore_hist.chrom{1}.png".format(
                self.sample, self.chrom_index)
        self.mapping.plot_hist_zscore(filename=self.directory + os.sep + 
                self.jinja["bp_plot"])

        # get region of interest
        roi = self.mapping.get_roi(self.features)

        # Low threshold case
        low_roi = roi.get_low_roi()
        low_cov_paragraph = ("Regions with a z-score lower than {0:.2f} and at "
            "least one base with a z-score lower than {1:.2f} are detected as "
            "low coverage region. Thus, there are {2} low coverage regions")
        self.jinja["lc_paragraph"] = low_cov_paragraph.format(
                    self.mapping.thresholds.low2,
                    self.mapping.thresholds.low, len(low_roi))

        # Save information relatd to the low ROIs
        html = self.htmltable(low_roi, "low_coverage", bgcolors=["size"])

        # Create a link for each row on the start position to jump directly to a
        # sub mapping page.
        step = self._get_step()
        def get_link(pos):
            start = divmod(pos, step)[0]
            name = "{0}_mapping_{1}".format(self.sample, start*step)
            stop = (start + 1)  * step
            if stop > len(self.mapping):
                stop = len(self.mapping)
            name = "{0}_{1}".format(name, stop)
            link = "submapping/{0}.chrom{1}.html".format(name, self.chrom_index)
            formatter = '<a target="_blank" href={0}>{1}</a>'
            return formatter.format(link, pos)
        html.to_csv()
        html.to_json()
        html.df["start"] = html.df["start"].apply(lambda x: get_link(x))
        self.jinja[html.name] = html.to_html()

        # Save information relatd to the high ROIs
        high_roi = roi.get_high_roi()
        high_cov_paragraph = ("Regions with a z-score higher than {0:.2f} and at "
            "least one base with a z-score higher than {1:.2f} are detected as "
            "high coverage region. Thus, there are {2} high coverage regions")
        self.jinja['hc_paragraph'] = high_cov_paragraph.format(
                self.mapping.thresholds.high2,
                self.mapping.thresholds.high, len(high_roi))

        html = self.htmltable(high_roi, "high_coverage", bgcolors=["size"])
        html.to_csv()
        html.to_json()
        html.df["start"] = html.df["start"].apply(lambda x: get_link(x))
        self.jinja[html.name] = html.to_html()


        # Sub mapping with javascript
        df = self._generate_submapping(high_roi, low_roi)
        html = HTMLTable(df)
        self.jinja['list_submapping'] = html.to_html(index=False) 
