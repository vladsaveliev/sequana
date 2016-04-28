# Import -----------------------------------------------------------------------

import os
import pandas as pd
from reports import HTMLTable
from sequana.report_main import BaseReport

# Class ------------------------------------------------------------------------

class VCFReport(BaseReport):
    """
    """
    def __init__(self, directory="report", **kargs):
        super(VCFReport, self).__init__(jinja_filename="vcf/index.html",
                directory=directory, output_filename="vcf_filter.html", **kargs)
        self.jinja['title'] = "VCF Report"

    def set_data(self, data):
        self.vcf = data

    def parse(self):
        csv_path = self.directory + os.sep + "vcf_filter.csv"
        self.vcf.to_csv(csv_path, index=False, header=True)
        self.jinja["csv_link"] = "vcf_filter.csv"

        html = HTMLTable(self.vcf)
        self.jinja['vcf_filter'] = html.to_html(index=False)
