"""


"""
from .report_main import BaseReport

# a utility from external reports package
from reports import HTMLTable

import pandas as pd


class AdapterRemovalReport(BaseReport):
    """A parent child to create report related to AdapterRemoval

    This can be used to design a new class for dedicated outputs from
    other tools such as for AlienTrimmer, CutAdapt, Skewer and so on.

    This class defines a minimal set of information to be provided

    """
    def __init__(self, jinja_template="adapter_removal",
            output_filename="adapter_removal.html", directory="report",
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
        super(AdapterRemovalReport, self).__init__(jinja_template, 
            output_filename, directory, **kargs)

        # Here, we defined default values from what is expected from the Jinja
        # template in share/templates/adapter_removal

        # generic: is it paired on single ?
        self.jinja['mode'] = 'unknown'

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
        # Create a Table with stats
        df = pd.DataFrame()
        df = pd.DataFrame({'Number of reads': [], 'percent': []})
        df.ix['Total reads'] = [
                    self.jinja['total_reads'],
                    '(100%)']
        df.ix['Too short'] = [
                    self.jinja['reads_too_short'],
                    self.jinja['reads_too_short_percent']]
        df.ix['Kept reads'] = [
                    self.jinja['reads_kept'],
                    self.jinja['reads_kept_percent']]

        h = HTMLTable(df)
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
        h = HTMLTable(df)
        html = h.to_html(index=True)
        self.jinja['adapters'] = html

        super(AdapterRemovalReport, self).create_report(onweb=onweb)















