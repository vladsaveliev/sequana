import easydev
import os


from .report_main import BaseReport

# a utility from external reports package
from reports import HTMLTable

import pandas as pd


def _get_template_path(name):
    # Is it a local directory ?
    if os.path.exists(name):
          return name
    else:
          template_path = easydev.get_shared_directory_path("sequana")
          template_path += os.sep + "templates"  + os.sep + name
          return template_path
  


class FixReport(BaseReport):
    """


    """
    def __init__(self, jinja_template="fix_contaminant",
            output_filename="fix.html", directory="report",
            overwrite=False, **kargs):
        """

        :param jinja_template: name of a directory (either local) or
            from sequana/share/templates where JINJA files are available. A file
            named index.html is required but may be renamed (with
            **output_filename** parameter).
        :param output_filename: name of the final HTML file.
        :param directory: name of the output directory (defaults to report)

        Parameters accepted by :class:`reports.Report` are also accepted.

        """
        super(FixReport, self).__init__(jinja_template, output_filename,
              directory, **kargs)

        self.title = "Fix Report Summary"
        self.jinja['title'] = "Fix Report Summary"

        self.jinja['mode'] = "Paired-end"
        self.mode = "pe"

        self.input_filename = "fix_stats.json"

    def parse(self):
        import json
        data = json.load(open(self.input_filename, "r"))

        for key, value in data.items():
            self.jinja[key] = value

        x = data['R1_mapped']
        y = data['R1_unmapped']

        contamination = x / float(x+y) * 100

        from easydev import precision
        self.jinja['contamination'] = precision(contamination, 3)
        

        import pandas as pd
        df = pd.DataFrame({
            'R1': [data['R1_mapped'], data['R1_unmapped']],
            'R2': [data['R2_mapped'], data['R2_unmapped']]})
        df.index = ['mapped', 'unmapped']
        
        
        h = HTMLTable(df)
        html = h.to_html(index=True)
        self.jinja['stats'] = html








