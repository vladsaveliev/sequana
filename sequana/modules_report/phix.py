# coding: utf-8
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
"""Module to write coverage report"""
import os
import glob
import io

from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils import config

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana import logger

from sequana.utils.datatables_js import DataTable


class PhixModule(SequanaBaseModule):
    """ Write HTML report of fastq stats analysis."""
    def __init__(self, input_directory, output_filename=None, tag_R1="_R1_"):
        """
        :param input_directory: where to find the json and boxplot image. The
            path where to find the data does not matter since the JSON and PNG
            will be embedded.
        :param path_to_fastqc: This must be provided by the user. This is the
            directory where will be found the original FastQC reports. This can
            be infered but is prone to error so for now, we must provide this
            argument.
        :param output_filename: if not provided, the HTML is not created.

        ::

            from sequana.modules_report.fastq_stats import FastQStatsModule
            ff = PhixModule("./SAMPLE/", "test.html")

        Expect to find bwa_mem_stats/bwa_mem_stats.json and phix_stats/*json

        where *json is a pattern SAMPLE_R1_.mapped.json or
        SAMPLE_R1_.unmapped.json. Same with _R2_

        """
        super().__init__()
        self.directory = input_directory
        self.tag_R1 = tag_R1
        self.filename1 = self.directory + "/bwa_bam_to_fastq/bwa_mem_stats.json"
        self.phix_directory = "fastq_stats_phix"
        self.create_report_content()
        if output_filename:
            self.create_html(output_filename)

    def create_report_content(self):
        """ Generate the sections list to fill the HTML report.
        """
        self.sections = list()
        self.add_stats()

    def _get_files(self, pattern):
        filenames = glob.glob(os.sep.join([self.directory, self.phix_directory, 
                                           pattern]))
        if len(filenames) == 4:
            mode = "pe"
        elif len(filenames) == 2:
            mode = "se"
        elif len(filenames) == 0:
            return
        else:
            logger.warning("PhixModule: more than 4 files "
                           "matched the pattern %s" % pattern)
            return
        return filenames, mode

    def _get_summary(self):
        from sequana.tools import StatsBAM2Mapped
        data = StatsBAM2Mapped(self.filename1).data
        return data

    def _get_html_summary_section(self):
        from easydev import precision
        data = self._get_summary()
        html = "Reads with Phix: %s %%<br>" % precision(data['contamination'], 3)
        html += "Unpaired: %s <br>" % data['unpaired']
        html += "duplicated: %s <br>" % data['duplicated']
        return html

    def _get_stats(self):
        filenames, mode = self._get_files("*.json")
        cols = ["A", "C", "G", "T", "GC content", "N", "average read length", 
                "mean quality", "n_reads", "total bases"]
        N = len(filenames)
        import numpy as np
        df = pd.DataFrame(np.zeros((N, 10)), columns=cols)

        indices = []
        for i, filename in enumerate(filenames):
            if self.tag_R1 in filename:
                index = "R1"
            else:
                index = "R2"
            if "unmapped" in filename:
                index+=".unmapped"
            else:
                index+=".mapped"
            indices.append(index)
            try:
                subdf = pd.read_json(filename)
                df.ix[i] = subdf.ix[0]
            except:
                pass
        df.index = indices
        return df

    def _get_html_stats_section(self):
        df = self._get_stats()
        datatable = DataTable(df, "phix_stats", index=True)
        datatable.datatable.datatable_options = {
            'scrollX': '300px',
            'pageLength': 15,
            'scrollCollapse': 'true',
            'dom': 'tpB',
            "paging": "false",
            'buttons': ['copy', 'csv']}
        js = datatable.create_javascript_function()
        html_tab = datatable.create_datatable(float_format='%.3g')
        html = """<p>We mapped the raw reads on a reference (see config file). The reads mapped are removed and the unmapped reads are kept. Here below are some statistics about the mapped and unmapped reads.
</p><p>
The A, C, G, T, N rows report the percentage of each bases in the overall
sequences</p>""" 
        html += "{} {}".format(html_tab, js)
        return html

    def _get_html(self):
        html1 = self._get_html_summary_section()
        html2 = self._get_html_stats_section()
        return html1 + html2

    def add_stats(self):
        self.sections.append({
            "name": "Stats inputs",
            "anchor": "phix_stats",
            "content": self._get_html()
        })

