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

from sequana.lazy import pandas as pd

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

    def get_mean_quality_samples(self):
        this = json.loads(self.data["sample_stats__samples_json"])
        return this["mean quality"]['0']

    def get_nreads_raw(self):
        this = json.loads(self.data["sample_stats__samples_json"])
        return this["n_reads"]['0']

    def get_average_read_length(self):
        this = json.loads(self.data["sample_stats__samples_json"])
        return this["average read length"]['0']

    def get_trimming_percent(self):
        d = self.get_cutadapt_stats()
        try:
            trimming = d["percent"]["Pairs too short"]
        except:
            trimming = d["percent"]["Reads too short"]
        return float(trimming.replace("%", "").replace("(","").replace(")","").strip())

    def get_adapters_percent(self):
        d = self.get_cutadapt_stats()
        try:
            read1 = d["percent"]["Read1 with adapters"]
            # e.g. (8.3%)
            read1 = float(read1.strip("%)").strip("("))
            try:
                read2 = d["percent"]["Read2 with adapters"]
                read2 = float(read2.strip("%)").strip("("))
                read1 = (read1 + read2 ) /2. # FIXME crude approximation
            except:
                pass # single-end
        except:
            read1 = d["percent"]["Reads with adapters"]
        return read1

    def get_read1_with_adapters_percent(self):
        return self._get_read_with_adapters_percent("1")
    def get_read2_with_adapters_percent(self):
        return self._get_read_with_adapters_percent("2")
    def _get_read_with_adapters_percent(self, tag):
        d = self.get_cutadapt_stats()
        trimming = d["percent"]["Read%s with adapters" % tag]
        trimming = trimming.strip()
        for this in [",", "(", ")", "%"]:
            trimming = trimming.replace(this, "")
        trimming = float(trimming)
        return trimming


class SequanaMultipleSummary(BaseReport):
    """
    Used by the pipelines to create a summary based on the content of the
    directory. Also used by the standalone application, in which case
    config and pipeline files are not required.

    For developers:

    1. In class Summary:

    2 In class SequanaMultipleSummary:

        try:self.populate_gc()
        except:pass
        def populate_gc():
            do something
        def get_gc(self):
            return [x.get_gc() for x in self.summaries]

    3: update the jinja file report_multiple_summay

    """
    def __init__(self, pattern="**/summary.json", verbose=True, **kargs):

        super(SequanaMultipleSummary, self).__init__(
            jinja_filename="multi_summary.html",
            directory=".",
            output_filename="multi_summary.html",
            **kargs)

        print("Sequana Summary is still a tool in progress and have been " +
              "  tested with the quality_control pipeline only for now.")
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

        self.jinja['links'] = [{'href': this.replace(".json", ".html"),
                                'caption': this.split("/",1)[0]}
                               for this in self.filenames]

        self.jinja['canvas'] = '<script type="text/javascript" src="js/canvasjs.min.js"></script>'
        self.jinja['canvas'] += """<script type="text/javascript">
            window.onload = function () {"""

        self.jinja['sections'] = []

        # The order does not matter here, everything is done in JINJA
        try:self.populate_nreads_raw()
        except Exception as err:
            print(err)

        try:self.populate_phix()
        except Exception as err:
            print(err)

        try:self.populate_gc_samples()
        except Exception as err:
            print(err)

        try: self.populate_trimming()
        except Exception as err:
            print(err)

        try:self.populate_mean_quality()
        except Exception as err:
            print(err)

        try:self.populate_adapters()
        except Exception as err:
            print(err)


        self.jinja['canvas'] += """
    function onClick(e){
        window.open(e.dataPoint.url)
    }
}</script>"""

    def _get_div(self, name, title):
        div = "<h2>%s</h2>" %  title
        div += '<div id="chartContainer' + name
        div += '" style="height: 300px; width: 90%;"></div></hr>'
        return div

    def _get_df(self, method_name):
        data = getattr(self, method_name)()
        df = pd.DataFrame({
            "name": self.get_projects(),
            "value": data,
            "url": self.get_urls()})
        return df

    def populate_adapters(self):
        title = "Adapters content"
        df = self._get_df("get_adapters_percent")
        cb = CanvasBar(df, "Adapters content", "adapters", xlabel="Percentage")
        self.jinja['canvas'] += cb.to_html()
        self.jinja['sections'].append(self._get_div("adapters", title))

    def populate_nreads_raw(self):
        df = self._get_df("get_nreads_raw")
        cb = CanvasBar(df, "Number of reads (raw data)", "nreads_raw",
                    xlabel="Number of reads")
        self.jinja['canvas'] += cb.to_html()
        self.jinja['sections'].append(self._get_div("nreads_raw",
                                                    "Number of reads"))

    def populate_mean_quality(self):
        title = "Mean quality"
        df = self._get_df("get_mean_quality_samples")
        cb = CanvasBar(df, title, "mean_quality", xlabel="mean quality")
        self.jinja['canvas'] += cb.to_html(options={'maxrange':40})
        self.jinja['sections'].append(self._get_div("mean_quality", title))

    def populate_gc_samples(self):
        title = "GC content"
        df = self._get_df("get_gc_content_samples")
        cb = CanvasBar(df, title, "populate_gc_samples", xlabel="Percentage")
        self.jinja['canvas'] += cb.to_html(options={"maxrange":100})
        self.jinja['sections'].append(self._get_div("populate_gc_samples",title))

    def populate_phix(self):
        title = "Phix content"
        df = self._get_df("get_phix_percent")
        cb = CanvasBar(df, title, "phix", xlabel="Percentage")
        self.jinja['canvas'] += cb.to_html()
        self.jinja['sections'].append(self._get_div("phix", title))

    def populate_trimming(self):
        title = "Trimming"
        df = self._get_df("get_trimming_percent")
        cb = CanvasBar(df, title, "trimming", xlabel="Percentage")
        self.jinja['canvas'] += cb.to_html()
        self.jinja['sections'].append(self._get_div("trimming", title))

    def get_projects(self):
        return [x.data['project'] for x in self.summaries]

    def get_urls(self):
        return [x.replace("summary.json", "summary.html") for x in self.get_unique_names()]

    def get_unique_names(self):
        """reduce the filenames length removing the common suffix
        remove also the filename itself.
        """
        projects = self.get_projects()
        if len(projects) == set(projects):
            return projects
        else:
            return self.filenames

    ##########################################################################
    # The retrieval of all specific data for all summaries

    def get_adapters_percent(self):
        return [x.get_adapters_percent() for x in self.summaries]

    def get_nreads_raw(self):
        return [x.get_nreads_raw() for x in self.summaries]

    def get_phix_percent(self):
        return [x.get_phix_percent() for x in self.summaries]

    def get_gc_content_samples(self):
        return [x.get_fastq_stats_samples()['GC content']["0"] for x in self.summaries]

    def get_trimming_percent(self):
        return [x.get_trimming_percent() for x in self.summaries]

    def get_mean_quality_samples(self):
        return [x.get_mean_quality_samples() for x in self.summaries]

    def parse(self):
        pass

