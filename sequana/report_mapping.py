# Import -----------------------------------------------------------------------

import os
import pandas as pd
from reports import HTMLTable
from sequana.report_main import BaseReport
from sequana.report_submapping import SubMappingReport

# Class ------------------------------------------------------------------------
class MappingReport(BaseReport):
    """
    """
    def __init__(self, low_threshold=-3, high_threshold=3, 
            directory="report", **kargs):
        super(MappingReport, self).__init__(
                jinja_filename="mapping/index.html",
                directory=directory,
                output_filename="report_mapping.html", **kargs)
        self.jinja['title'] = "Mapping Report"
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
            name = "mapping_{0}".format(i)
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
        self.jinja['main_link'] = 'index.html'

        self.mapping.plot_coverage(filename=self.directory + os.sep + 
                "coverage.png", low_threshold=self.low_t, 
                high_threshold=self.high_t)
        
        self.jinja["mu"] = "{:.2f}".format(self.mapping.best_gaussian["mu"])
        self.jinja["sigma"] = "{:.2f}".format(
                self.mapping.best_gaussian["sigma"])
        self.mapping.plot_hist(filename=self.directory + os.sep +
                "zscore_hist.png")

        low_cov_df = self.mapping.get_low_coverage(self.low_t / 2)
        merge_low_cov = low_cov_df.merge_region(self.low_t)
        self.jinja['low_cov_threshold'] = self.low_t
        self.jinja["low_cov_threshold_2"] = "{0:.2f}".format(
                float(self.low_t) / 2)
        self.jinja['nb_low_region'] = len(merge_low_cov)
        html = HTMLTable(merge_low_cov)
        html.add_bgcolor("size")
        self.jinja['low_coverage'] = html.to_html(index=False)
        
        high_cov_df = self.mapping.get_high_coverage(self.high_t / 2)
        merge_high_cov = high_cov_df.merge_region(self.high_t)
        self.jinja['high_cov_threshold'] = self.high_t
        self.jinja["high_cov_threshold_2"] = "{0:.2f}".format(
                float(self.high_t) / 2)
        self.jinja['nb_high_region'] = len(merge_high_cov)
        html = HTMLTable(merge_high_cov)
        html.add_bgcolor("size")
        self.jinja['high_coverage'] = html.to_html(index=False)

        df = self._generate_submapping(merge_low_cov, merge_high_cov)
        html = HTMLTable(df)
        self.jinja['list_submapping'] = html.to_html(index=False) 
