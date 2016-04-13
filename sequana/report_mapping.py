# Import -----------------------------------------------------------------------

import os
from reports import HTMLTable
from sequana.report_main import BaseReport

# Class ------------------------------------------------------------------------
class MappingReport(BaseReport):
    """
    """
    def __init__(self, low_threshold=-3, high_threshold=3, **kargs):
        super(MappingReport, self).__init__(
                jinja_filename="mapping/index.html",
                directory="report",
                output_filename="mapping.html", **kargs)
        self.jinja['title'] = "Mapping Report"
        self.low_t = low_threshold
        self.high_t = high_threshold

    def set_data(self, data):
        self.mapping = data

    def parse(self):
        self.mapping.plot_coverage(filename=self.directory + os.sep + 
                                            "coverage.png")
        self.mapping.plot_hist(filename=self.directory + os.sep + 
                                            "zscore_hist.png")
        low_cov_df = self.mapping.get_low_coverage(self.low_t)
        merge_low_cov = self.mapping.merge_region(low_cov_df)
        html = HTMLTable(merge_low_cov)
        html.add_bgcolor("size")
        self.jinja['low_coverage'] = html.to_html(index=False)
        
        high_cov_df = self.mapping.get_high_coverage(self.high_t)
        merge_high_cov = self.mapping.merge_region(high_cov_df)
        html = HTMLTable(merge_high_cov)
        html.add_bgcolor("size")
        self.jinja['high_coverage'] = html.to_html(index=False)
