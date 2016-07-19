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

import pandas as pd
import pysam

from sequana.reports.report_main import BaseReport
from sequana.bamtools import SAMFlags

from reports import HTMLTable


class MappingReport(BaseReport):
    """Report dedicated to Mapping

    ::

        from sequana import bedtools, sequana_data
        from sequana.reports.report_mapping import MappingReport
        mydata = bedtools.GenomeCov(sequana_data("test_bedcov.bed"))

        r = MappingReport()
        r.set_data(mydata)
        r.create_report() 

    """
    def __init__(self, directory="report", project="",
            output_filename="report_mapping.html", **kargs):
        super(MappingReport, self).__init__(
                jinja_filename="mapping/index.html", directory=directory,
                output_filename=output_filename, **kargs)
        self.project = project
        self.jinja['title'] = "Mapping Report of {0}".format(project)

    def set_data(self, chrom_list, bam=None):
        self.chrom_list = chrom_list
        self.bam = bam

    def parse(self):        
        self.jinja['main_link'] = 'index.html'
        self.jinja['bam_is_present'] = True

        try:
            self.jinja['alignment_count'] = len(self.bam)

            # first, we store the flags
            df = self.bam.get_flags_as_df().sum()
            df = df.to_frame()
            df.columns = ['counter']
            sf = SAMFlags()
            df['meaning'] = sf.get_meaning()
            df = df[['meaning', 'counter']]
            html = HTMLTable(df).to_html(index=True)
            self.jinja['flags_table'] = html

            # create the bar plot with flags
            self.bam.plot_bar_flags(logy=True, filename=self.directory + 
                    os.sep + "bar_flags_logy.png")
            self.bam.plot_bar_flags(logy=False, filename=self.directory + 
                    os.sep + "bar_flags.png")
            self.bam.plot_bar_mapq(filename=self.directory + os.sep + 
                    "bar_mapq.png")
        except TypeError:
            self.jinja['bam_is_present'] = False

        df = pd.DataFrame()
        formatter = '<a target="_blank" alt={0} href="{1}">{0}</a>'
        for chrom in self.chrom_list:
            link = self.project + "_" + chrom.chrom_name + "_mapping.html"
            df = df.append({"chromosome": formatter.format(
                chrom.chrom_name, link), "size": "{0:,}".format(len(chrom))}, 
                ignore_index=True)
        html = HTMLTable(df)
        self.jinja['list_chromosome'] = html.to_html(index=False)
