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
  


class FastQCReport(BaseReport):
    """


    """
    def __init__(self, output_filename="fastqc.html", directory="report", **kargs):
        """

        :param jinja_template: name of a directory (either local) or
            from sequana/share/templates where JINJA files are available. A file
            named index.html is required but may be renamed (with
            **output_filename** parameter).
        :param output_filename: name of the final HTML file.
        :param directory: name of the output directory (defaults to report)

        Parameters accepted by :class:`reports.Report` are also accepted.

        """
        super(FastQCReport, self).__init__(
                jinja_filename="fastqc/index.html", 
                directory=directory, 
                output_filename=output_filename, **kargs)


        self.jinja['title'] = "FastQC Report Summary"
        self.jinja['mode'] = "Paired-end"
        self.mode = "pe"
        self.input_filename = "fastqc.json"

    def parse(self):
        # Where to find the data ?
        # Assuming in fastqc/ directory, we figure out the HTML files and
        # create a table accordingly.

        import glob
        links = glob.glob("fastqc/*html")

        names = [filename.split('/')[1].split('.html')[0] for filename in links]

        df = pd.DataFrame({"names": names})
        df.sort_values(by='names')

        formatter = '<a target="_blank" alt={1} href="../{0}">{1}</a>'
        df["names"] = [formatter.format(link,name) for link,name in zip(links, names)]

        h = HTMLTable(df)
        html = h.to_html(index=True)
        self.jinja['list_fastqc'] = html








