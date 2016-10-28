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
import glob
import json

from sequana.reporting.report_main import BaseReport
from easydev import DevTools
from sequana.resources.canvas.bar import CanvasBar

import pandas as pd

__all__ = ['SequanaMultipleSummary']


class ReadSummary(object):
    def __init__(self, filename):
        self.filename = filename
        self.data = json.load(open(self.filename, "r"))

    def get_phix_percent(self):
        return self.data['phix_section_json']['contamination']

    def get_cutadapt_stats(self):
        return json.loads(self.data["cutadapt_stats1_json"])

    def get_fastq_stats_samples(self):
        return json.loads(self.data["sample_stats__samples_json"])

    def get_mean_quality__samples(self):
        return json.loads(self.data["sample_stats__samples_json"])["mean quality"]['0']

    def get_trimming_percent(self):
        d = self.get_cutadapt_stats()
        trimming = d["percent"]["Pairs too short"]
        return float(trimming.replace("%", "").replace("(","").replace(")","").strip())

    def get_read1_with_adapters_percent(self):
        d = self.get_cutadapt_stats()
        trimming = d["percent"]["Read1 with adapters"]
        return int(trimming.replace(",", "").strip())

    def get_read2_with_adapters_percent(self):
        d = self.get_cutadapt_stats()
        trimming = d["percent"]["Read2 with adapters"]
        return int(trimming.replace(",", "").strip())


class SequanaMultipleSummary(BaseReport):
    """
    Used by the pipelines to create a summary based on the content of the
    directory. Also used by the standalone application, in which case
    config and pipeline files are not required.


    """
    def __init__(self, pattern="**/summary.json", verbose=True, **kargs):

        super(SequanaMultipleSummary, self).__init__(
            jinja_filename="multi_summary.html",
            directory=".",
            output_filename="multi_summary.html",
            **kargs)

        print("Sequana Summary is still a tool in progress and would only work with the output of the quality_control .")
        self.verbose = verbose
        workdir = "."
        self.jinja['title'] = "Sequana multiple summary" 
        self.env.loader.searchpath.extend([workdir])
        self.devtools = DevTools()

        self.filenames = list(glob.iglob(pattern, recursive=True))
        self.summaries = [ReadSummary(filename) for filename in self.filenames]

        if self.verbose:
            print("Found %s projects/samples/ directories" % len(self.summaries))

        # The base has a navigation, that we do not want
        self.jinja['nav_off'] = 'True'

        self.jinja['n_samples'] = len(self.summaries)

        self.jinja['canvas'] = '<script type="text/javascript" src="js/canvasjs.min.js"></script>'
        self.jinja['canvas'] += """<script type="text/javascript">
            window.onload = function () {"""

        try:self.populate_phix()
        except:pass
        try:self.populate_gc_samples()
        except:pass
        try: self.populate_trimming()
        except:pass
        try:self.populate_mean_quality()
        except:pass

        self.jinja['canvas'] += "}</script>"

    def populate_mean_quality(self):
        datadict = dict([(k, v) for k,v in zip(self.get_unique_names(), 
            self.get_mean_quality__samples())])
        cb = CanvasBar(datadict, "Mean quality", "_mean_quality", xlabel="mean quality")
        self.jinja['canvas'] += cb.to_html()
        self.jinja['mean_quality'] =  '<div id="chartContainer_mean_quality" style="height: 300px; width: 100%%;">'

    def populate_gc_samples(self):
        datadict = dict([(k, v) for k,v in zip(self.get_unique_names(), 
            self.get_gc_content_samples())])
        cb = CanvasBar(datadict, "GC content", "populate_gc_samples", xlabel="Percentage")
        self.jinja['canvas'] += cb.to_html()
        self.jinja['gc'] =  '<div id="chartContainerpopulate_gc_samples" style="height: 300px; width: 100%%;">'

    def populate_phix(self):
        # Phix content
        datadict = dict([(k, v) for k,v in zip(self.get_unique_names(), self.get_phix_percent())])
        from sequana.resources.canvas.bar import CanvasBar
        # dictionary is not sorted...
        cb = CanvasBar(datadict, "Phix content", "phix", xlabel="Percentage")
        self.jinja['canvas'] += cb.to_html() 
        self.jinja['phix'] =  '<div id="chartContainerphix" style="height: 300px; width: 100%%;">'

    def populate_trimming(self):
        # Phix content
        datadict = dict([(k, v) for k,v in zip(self.get_unique_names(),
            self.get_trimming_percent())])
        from sequana.resources.canvas.bar import CanvasBar
        # dictionary is not sorted...
        cb = CanvasBar(datadict, "Trimming", "trimming", xlabel="Percentage")
        self.jinja['canvas'] += cb.to_html() 
        self.jinja['trimming'] =  '<div id="chartContainertrimming" style="height: 300px; width: 100%%;">'

    def get_cutadapt_stats1(self):
        return [x.get_cutadapt_stats() for x in self.summaries]

    def get_phix_percent(self):
        return [x.get_phix_percent() for x in self.summaries]

    def get_projects(self):
        return [x.data['project'] for x in self.summaries]

    def get_unique_names(self):
        """reduce the filenames length removing the common suffix
        remove also the filename itself.
        """
        projects = self.get_projects()
        if len(projects) == set(projects):
            return projects
        else:
            return self.filenames

    def get_gc_content_samples(self):
        return [x.get_fastq_stats_samples()['GC content']["0"] for x in self.summaries]

    def get_trimming_percent(self):
        return [x.get_trimming_percent() for x in self.summaries]

    def get_mean_quality__samples(self):
        return [x.get_mean_quality__samples() for x in self.summaries]

    def parse(self):
        pass

