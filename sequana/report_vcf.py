# Import -----------------------------------------------------------------------

import os
import shutil
import pandas as pd
from reports import HTMLTable
from sequana.report_main import BaseReport

# Class ------------------------------------------------------------------------

class VCFReport(BaseReport):
    """
    """
    def __init__(self, csv_file, directory="report", **kargs):
        super(VCFReport, self).__init__(jinja_filename="vcf/index.html",
                directory=directory, output_filename="report_vcf.html", **kargs)
        self.jinja['title'] = "VCF Report"
        self.csv = csv_file

    def set_data(self, data):
        self.vcf_record = data

    def parse(self):
        df = self.vcf_record.vcf_to_csv(self.csv)
        self.jinja["csv_link"] = self.csv.split("/")[-1]

        shutil.copy(self.vcf_record.filename, self.directory)
        self.jinja["vcf_link"] = self.vcf_record.filename.split("/")[-1] 

        html = HTMLTable(df)
        self.jinja['vcf_filter'] = html.to_html(index=False)
