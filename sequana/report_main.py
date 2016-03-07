import os
import easydev

from reports import Report
from sequana import version


def _get_template_path(name):
    # Is it a local directory ?
    if os.path.exists(name):
        return name
    else:
        template_path = easydev.get_shared_directory_path("sequana")
        template_path += os.sep + "templates"  + os.sep + name
        return template_path


class BaseReport(Report):
    """A Parent child for all reports created in Sequana



    """
    def __init__(self, jinja_template, output_filename,
                 directory="report", **kargs):
        """.. rubric:: Constructor

        :param jinja_template: name of a directory (either local) or
            from sequana/share/templates where JINJA files are available. A file
            named index.html is required but may be renamed (with
            **output_filename** parameter).
        :param output_filename: name of the final HTML file.
        :param directory: name of the output directory (defaults to report)

        Parameters accepted by :class:`reports.Report` are also accepted.

        """
        # finds automatically the local directory or a directory to be found in
        # sequana distribution
        template_path = _get_template_path(jinja_template)
        super(BaseReport, self).__init__(template_path=template_path,
            directory=directory, **kargs)

        # This redefines the default name of the output (index.html) otherwise,
        # several reports will overwrite the default index.html. One must make
        # sure that adapter_removal does not exists
        self.filename = output_filename

        # That is where we will store all data to be used by the Jinja templates
        # in Report, everything is saved in jinja, but just to not get confused,
        # we will use data attribute for now
        self.data = {}

        # Here, we defined default values from what is expected from the Jinja
        # template in share/templates/adapter_removal
        self.data['sequana_version'] = version

        # Common information to be filled (possibly)
        self.data['command'] = "unset"

    def parse(self):
        """populate the :attr:`data` attribute used by the JINJA templates

        """
        raise NotImplementedError

    def create_report(self, onweb=False):
        try:
            self.parse()
        except Exception as err:
            print(err)
            print("Parsing of %s failed" % self.input_filename)

        super(BaseReport, self).create_report(onweb=onweb)


class SequanaReport(BaseReport):
    def __init__(self, jinja_template="main", output_filename="index.html", 
        directory="report", **kargs):

        super(SequanaReport, self).__init__(jinja_template, output_filename,
            directory, **kargs)

        self.title = "Sequana Report"

        #self.data['snakefile'] = "undefined"

    def read_snakefile(self, filename):
        """

        :param str filename:
        """
        with open(filename, "r") as fin:
            self.jinja['snakefile'] = fin.read()

    def parse(self, snakefile="Snakefile"):
        self.read_snakefile(snakefile)
        










