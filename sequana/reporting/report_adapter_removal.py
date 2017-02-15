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
"""


"""
import os
from sequana.reporting.report_main import BaseReport

# a utility from external reports package
from sequana.lazy import reports

from sequana.lazy import pandas as pd
from sequana.lazy import pylab


class AdapterRemovalReport(BaseReport):
    """A parent child to create report related to AdapterRemoval

    This can be used to design a new class for dedicated outputs from
    other tools such as for AlienTrimmer, CutAdapt, Skewer and so on.

    This class defines a minimal set of information to be provided

    """
    def __init__(self,
            output_filename="adapter_removal.html", 
            directory="report",
            overwrite=False, **kargs):
        """

        :param jinja_template: name of a directory (either local) or
            from sequana/share/templates where JINJA files are available. A file
            named index.html is required but may be renamed (with
            **output_filename** parameter).
        :param output_filename: name of the final HTML file.
        :param directory: name of the output directory (defaults to report)

        Parameters accepted by :class:`sequana.reporting.Report` are also accepted.

        """
        super(AdapterRemovalReport, self).__init__(
            jinja_filename="adapter_removal/index.html",
            output_filename=output_filename,
            directory=directory, **kargs)

        # Here, we defined default values from what is expected from the Jinja
        # template in share/templates/adapter_removal

        # some stats
        self.jinja['total_reads'] = None
        self.jinja["total_reads_percent"] = None
        self.jinja["reads_too_short"] = None
        self.jinja["reads_too_short_percent"] = None
        self.jinja['reads_kept'] = None
        self.jinja['reads_kept_percent'] = None

        # adapters
        self.jinja['adapters'] = []

        # should be filled with a table with basic stats
        self.jinja['stats'] = ""

    def create_report(self, onweb=False):

        self.parse()
        df = pd.DataFrame()

        if self.mode == "pe":
            prefix = "paired_"
        else:
            prefix = ""

        df = pd.DataFrame({'Number of reads': [], 'percent': []})
        df.ix['Total paired reads'] = [
                    self.jinja['%stotal_reads' % prefix],
                    '(100%)']
        if self.mode == "pe":
            df.ix['Read1 with adapters'] = [
                    self.jinja['%sreads1_with_adapters' % prefix],
                    self.jinja['%sreads1_with_adapters_percent'% prefix]]
            df.ix['Read2 with adapters'] = [
                    self.jinja['%sreads2_with_adapters' % prefix],
                    self.jinja['%sreads2_with_adapters_percent'% prefix]]
        else:
            df.ix['Pairs with adapters'] = [
                    self.jinja['%sreads_with_adapters' % prefix],
                    self.jinja['%sreads_with_adapters_percent'% prefix]]
        df.ix['Pairs too short'] = [
                    self.jinja['%sreads_too_short' % prefix],
                    self.jinja['%sreads_too_short_percent'% prefix]]
        df.ix['Pairs kept'] = [
                    self.jinja['%sreads_kept' % prefix],
                    self.jinja['%sreads_kept_percent' % prefix]]
        if self.mode != "pe":
            df.index = [this.replace("paired", "").replace("Pairs", "Reads") for this in df.index]

        df.to_json(self.sample_name + "/cutadapt/cutadapt_stats1.json")

        h = reports.HTMLTable(df)
        html = h.to_html(index=True)
        self.jinja['stats'] = html

        # Create a Table with adapters
        df = pd.DataFrame()
        df = pd.DataFrame({'Length': [], 'Trimmed':[], 'Type':[], 'Sequence': [], })

        for count, adapter in enumerate(self.data['adapters']):
            name = adapter['name']
            info = adapter['info']
            df.ix[name] = [info['Length'], info['Trimmed'],
                info['Type'], info['Sequence']]
        df.columns = ['Length', 'Trimmed', 'Type', 'Sequence']

        df.to_json(self.sample_name + "/cutadapt/cutadapt_stats2.json")
        h = reports.HTMLTable(df)
        html = h.to_html(index=True)
        self.jinja['adapters'] = html


        # This would work for cutadapt only
        histograms = self.get_histogram_data()
        html = ""
        html += "<div>\n"
        from easydev import DevTools
        DevTools().mkdir(self.sample_name + "/cutadapt/images")
        for key in sorted(histograms.keys()):
            if len(histograms[key]) <= 1:
                continue
            histograms[key].plot(logy=True, lw=2, marker="o")
            pylab.title(name)
            name = key.replace(" ", "_")
            filename =  "%s/cutadapt/images/%s.png" % (self.sample_name,name)
            try:
                pylab.savefig(filename)
            except FileNotFoundError:
                print("Warning:: cutadapt report, image not created")
            pylab.grid(True)
            html += '<img src="cutadapt/images/%s.png" width="45%%"></img> ' % (name)
            #except:
            #    pass
        html += "</div>\n"

        self.jinja['cutadapt'] = html


        super(AdapterRemovalReport, self).create_report(onweb=onweb)















