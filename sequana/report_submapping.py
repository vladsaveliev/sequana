# Import -----------------------------------------------------------------------

import os
from reports import HTMLTable
from sequana.report_main import BaseReport

# Class ------------------------------------------------------------------------
class SubMappingReport(BaseReport):
    """
    """
    def __init__(self, start, stop, low_threshold=-3, high_threshold=3, 
            directory="report", output_filename="submapping.html", **kargs):
        super(SubMappingReport, self).__init__(
                jinja_filename="submapping/index.html",
                directory=directory,
                output_filename=output_filename, **kargs)
        self.jinja['title'] = "Mapping Report [{0},{1}]".format(start, stop)
        self.start = start
        self.stop = stop
        self.low_t = low_threshold
        self.high_t = high_threshold

    def set_data(self, data):
        self.mapping = data

    def parse(self):
        self.mapping.write_csv(self.directory + os.sep +
                "mapping_{0}_{1}.csv".format(self.start, self.stop), 
                start=self.start, stop=self.stop, header=False)
        self.jinja['input_df'] = "'mapping_{0}_{1}.csv'".format(self.start, 
                self.stop)

        low_cov_df = self.mapping.get_low_coverage(threshold=self.low_t, 
                start=self.start, stop=self.stop)
        merge_low_cov = self.mapping.merge_region(low_cov_df)
        self.jinja["low_cov_threshold"] = self.low_t
        self.jinja["nb_low_region"] = len(merge_low_cov)
        html = HTMLTable(merge_low_cov)
        html.add_bgcolor("size")
        self.jinja['low_coverage'] = html.to_html(index=False)
        
        high_cov_df = self.mapping.get_high_coverage(threshold=self.high_t, 
                start=self.start, stop=self.stop)
        merge_high_cov = self.mapping.merge_region(high_cov_df)
        self.jinja["high_cov_threshold"] = self.high_t
        self.jinja["nb_high_region"] = len(merge_high_cov)
        html = HTMLTable(merge_high_cov)
        html.add_bgcolor("size")
        self.jinja['high_coverage'] = html.to_html(index=False)
