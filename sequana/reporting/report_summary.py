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
from sequana.snaketools import SequanaConfig
from sequana.snaketools import Module
from sequana.snaketools import FastQFactory
from sequana import tools

from easydev import DevTools
from sequana.lazy import reports

from sequana.lazy import pandas as pd


__all__ = ['SequanaSummary']


class SequanaSummary(BaseReport):
    """
    Used by the pipelines to create a summary based on the content of the
    directory. Also used by the standalone application, in which case
    config and pipeline files are not required.


    Note that you may provide a manager to tell if paired or not.

    """
    def __init__(self,  sample, directory="report", output_filename="summary.html",
                    configfile="config.yaml", snakefile=None,
                    workdir=".", workflow=True, include_all=True,
                    manager=None, **kargs):

        super(SequanaSummary, self).__init__(
            jinja_filename="summary.html",
            directory=directory,
            output_filename=output_filename,
            **kargs)

        self.env.loader.searchpath.extend([workdir])

        self.workdir = workdir
        self.devtools = DevTools()
        self.sample = sample

        self.title = "Summary Report"
        self.jinja['title'] = "Summary report"

        self.manager = manager

        if self.manager.paired is True:
            self.jinja['type'] = "Paired-end"
        else:
            self.jinja['type'] = "Single-end"
        # ============================================ Add the config and pipeline files

        if configfile:
            self.config = SequanaConfig(configfile)
            self.jinja['project'] = sample

            try:
                pipeline_name = snakefile.split(".rules")[0]
                url = "http://sequana.readthedocs.io/en/master/pipelines.html#" + pipeline_name
                self.jinja["pipeline_name"] = '<a href="%s"> %s</a>' % (url,pipeline_name.title())
                self.jinja["pipeline_name"] += " -- <i>[%s]</i>" % Module(pipeline_name).overview
            except:
                # for the standalone apps
                pass

        # The base has a navigation, that we do not want
        self.jinja['nav_off'] = 'True'
        self.jinja['workflow'] = workflow

        if snakefile: self.read_snakefile(snakefile)

        # This is a string representation of the config file
        if configfile:
            try:self.read_configfile(configfile)
            except:pass

        # include whatever is relevant
        if include_all:
            try: self.include_kraken()
            except: pass

            try:self.include_phix()
            except:pass

            try: self.include_sample_stats()
            except:pass

            try:self.include_adapters_stats()
            except:pass

            self.include_details()
            self.include_input_links()
            self.include_output_links()

        # this is a dictionary usable within JINJA templates, which may have
        # been chnaged by methods above
        if configfile:
            self.jinja['cfg'] = self.config.config

        with open("%s/summary.json" % directory, "w") as fh:
            json.dump(self.jinja, fh, indent=4)

    def include_input_links(self):
        # Links to the datasets
        html = "<ul>"
        for fullpath in self.manager.samples[self.sample]:
            filename = os.path.basename(fullpath)
            html += '<li>Raw data: <a href="%s">%s</a></li>\n' % (fullpath, filename)
        html += "</ul>"
        self.jinja['dataset'] = html

    def include_output_links(self):
        html = "<ul>"
        filenames = glob.glob(self.directory+"/*fastq.gz")
        if len(filenames):
            ff = FastQFactory(filenames)
            for filename in ff.basenames:
                html += '<li>Download cleaned data: <a href="%s">%s</a></li>\n' % (filename,
                    filename)
            html += "</ul>"
            self.jinja['output'] = html

        # if cutadapt. If not TODO
        filenames = glob.glob("%s//fastq_stats_cutadapt/*boxplot.png" % self.directory)
        for filename in filenames:
            filename = filename.split("//", 1)[1].strip("/")
            if "R1" in filename:
                self.jinja['output_image_r1'] = filename
                self.jinja['output_image_r1_href'] = filename.replace(
                    "fastq_stats_cutadapt", "fastqc_cutadapt").replace(
                    "boxplot.png", "fastqc.html")
            elif "R2" in filename:
                self.jinja['output_image_r2'] = filename
                self.jinja['output_image_r2_href'] = filename.replace(
                    "fastq_stats_cutadapt", "fastqc_cutadapt").replace(
                    "boxplot.png", "fastqc.html")

    def include_details(self):
        self.jinja['snakemake_stats'] = "snakemake_stats.png"

    def include_adapters_stats(self):
        filename = self.directory + "/fastq_stats_cutadapt/temp.html"
        try:
            self.jinja['cutadapt_stats2'] = open(filename, "r").read()
        except:
            self.config.config['adapter_removal']['do'] = False

        if self.config.config["adapter_removal"]['do']:
            df = pd.read_json(self.directory + "/cutadapt/cutadapt_stats1.json")
            self.jinja["cutadapt_stats1_json"] = df.to_json()
            h = reports.HTMLTable(df)
            html = h.to_html(index=True)
            self.jinja['cutadapt_stats1'] = html

    def include_sample_stats(self):
        filename = self.directory + "/fastq_stats_samples/temp.html"
        self.jinja['sample_stats'] = open(filename, "r").read()

        filenames = glob.glob("%s//fastq_stats_samples/*boxplot.png" % self.directory)
        for filename in filenames:
            filename = filename.split("//", 1)[1].strip("/")
            if "R1" in filename:
                self.jinja['sample_image_r1'] = filename
                self.jinja['sample_image_r1_href'] = filename.replace(
                    "fastq_stats_samples", "fastqc_samples").replace(
                    "boxplot.png", "fastqc.html")
                newfilename  = self.directory + "/fastq_stats_samples/%s._R1_.json" % self.sample
                df = pd.read_json(newfilename)
                self.jinja["sample_stats__samples_json"] = df.to_json()

            elif "R2" in filename:
                self.jinja['sample_image_r2'] = filename
                self.jinja['sample_image_r2_href'] = filename.replace(
                    "fastq_stats_samples", "fastqc_samples").replace(
                    "boxplot.png", "fastqc.html")
                newfilename  = self.directory + "/fastq_stats_samples/%s._R2_.json" % self.sample
                df = pd.read_json(newfilename)
                self.jinja["sample_stats__samples_json_R2"] = df.to_json()

    def include_kraken(self):
        self.jinja['kraken_pie'] = "kraken/kraken.png"
        try:
            self.jinja['kraken_database'] = os.path.basename(
                self.config.config['kraken']['database'])
        except:
            self.jinja['kraken_database'] = "?"

        df = pd.read_csv(self.directory + "/kraken/kraken.csv")

        # Rounding and convert in string to avoid exp notation
        df['percentage']  = df['percentage'].apply(lambda x: str(round(x,4)))
        self.jinja['kraken_json'] = df.to_json()

        table = self.htmltable(df, tablename="kraken")
        if "ena" in table.df.columns:
            table.add_href('ena', url="http://www.ebi.ac.uk/ena/data/view/")
        table.name = "kraken/kraken"
        self.jinja['kraken_html_table'] = table.to_html(index=False)


    def include_phix(self):
        filename=self.directory + "/bwa_bam_to_fastq/bwa_mem_stats.json"
        if os.path.exists(filename):
            stats = tools.StatsBAM2Mapped(filename)
            self.jinja['phix_section'] = stats.to_html(with_stats=False)
            self.jinja['phix_section_json'] = stats.data
        else:
            print('Could not find phix information from file %s ' % filename)

        # include html pages with some stats
        filename = self.directory + "/fastq_stats_phix/temp.html"
        try:
            self.jinja['phix_stats'] = open(filename, "r").read()
        except:
            pass

    def parse(self):
        pass

