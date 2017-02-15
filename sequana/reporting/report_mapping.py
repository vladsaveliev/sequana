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
"""Report dedicated to Mapping"""
import os

from sequana.lazy import pandas as pd

from sequana.reporting.report_main import BaseReport
from sequana.lazy import reports


class MappingReport(BaseReport):
    """Report dedicated to Mapping

    ::

        from sequana import bedtools, sequana_data
        from sequana.reporting.report_mapping import MappingReport
        mydata = bedtools.GenomeCov(sequana_data("test_bedcov.bed"))

        r = MappingReport()
        r.set_data(mydata)
        r.create_report()

    """
    def __init__(self, directory="report", sample="",
                 output_filename="report_mapping.html", **kargs):
        super(MappingReport, self).__init__(
            jinja_filename="mapping/index.html", directory=directory,
            output_filename=output_filename, **kargs)
        self.sample = sample
        self.jinja['title'] = "Mapping Report of {0}".format(sample)

    def set_data(self, chrom_list, bam=None, quast=None):
        self.chrom_list = chrom_list
        self.bam = bam
        self.quast = quast

    def parse(self):
        self.jinja['main_link'] = 'index.html'

        if self.bam:
            self.jinja['bam_is_present'] = True

            # first, we store the flags
            df = self.bam.get_flags_as_df().sum()
            nb_mapped_read = len(self.bam) - (df.loc[4] + df.loc[256] +
                                              df.loc[2048])
            html = ("<ul>\n"
                    "<li>Reads mapped: {0}</li>\n"
                    "<li>Reads mapped in a proper pair: {1}</li>\n"
                    "<li>Reads unmapped: {2}</li>\n"
                    "<li>Reads with supplementary alignment (hard clipped): "
                    "{3}</li>\n"
                    "<li>Reads duplicated: {4}</li>\n"
                    "</ul>\n")
            self.jinja["summary_bam"] = html.format(nb_mapped_read, df.loc[2],
                                                    df.loc[4], df.loc[2048],
                                                    df.loc[1024])

            # create the bar plot with flags
            image_prefix = "images/" + self.sample
            self.jinja["bar_log_flag"] = image_prefix + "_bar_flags_logy.png"
            self.bam.plot_bar_flags(logy=True, filename=self.directory +
                                    os.sep + self.jinja["bar_log_flag"])
            self.jinja["bar_flag"] = image_prefix + "_bar_flags.png"
            self.bam.plot_bar_flags(logy=False, filename=self.directory +
                                    os.sep + self.jinja["bar_flag"])
            self.jinja["bar_mapq"] = image_prefix + "_bar_mapq.png"
            self.bam.plot_bar_mapq(filename=self.directory + os.sep +
                                   self.jinja["bar_mapq"])

        formatter = '<a target="_blank" alt={0} href="{1}">{0}</a>'
        if self.quast:
            self.jinja['title'] = "Denovo Report of {0}".format(self.sample)
            self.jinja['quast_is_present'] = True
            quast_link = formatter.format("here", self.quast)
            self.jinja['quast_link'] = ("The report provides by quast " +
                                        "is {0}.".format(quast_link))

        # create table with links to chromosome reports
        link = "{0}_mapping.chrom{1}.html"
        df = pd.DataFrame([[formatter.format(chrom.chrom_name,
                          link.format(self.sample, chrom.chrom_index)),
                          chrom.get_size(), chrom.get_mean_cov(),
                          chrom.get_var_coef()] for chrom in self.chrom_list],
                          columns=["chromosome", "size", "mean_coverage",
                          "coef_variation"])
        html = reports.HTMLTable(df)
        self.jinja['list_chromosome'] = html.to_html(index=False)
