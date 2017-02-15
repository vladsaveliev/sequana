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
"""Data structure to build reports in a consistent way"""
import os
import shutil


import easydev
sepjoin = os.sep.join

from sequana.lazy import reports
from sequana import version

import glob


class BaseReport(reports.Report):
    """A Parent child for all reports created in Sequana


    """
    def __init__(self, jinja_filename, directory="report",
                output_filename="test.html", init_report=True, **kargs):
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
        sequana_path = easydev.get_package_location('sequana')
        extra_css_path = sepjoin([sequana_path, "sequana", "resources", "css"])
        extra_js_path = sepjoin([sequana_path, "sequana", "resources", "js"])

        extra_css_list = glob.glob(extra_css_path + os.sep + "*css")
        extra_js_list = glob.glob(extra_js_path + os.sep + "*js")

        searchpath = sepjoin([sequana_path, "sequana", "resources", "jinja"])

        super(BaseReport, self).__init__(searchpath, filename=output_filename,
            template_filename=jinja_filename,
            directory=directory,
            extra_css_list=extra_css_list,
            extra_js_list=extra_js_list,
            init_report=init_report
            )

        # This redefines the default name of the output (index.html) otherwise,
        # several reports will overwrite the default index.html.
        self.filename = output_filename

        # That is where we will store all data to be used by the Jinja templates
        # in Report, everything is saved in jinja, but just to not get confused,
        # we will use data attribute for now
        #self.data = {}

        # Here, we defined default values from what is expected from the Jinja
        # template in share/templates/adapter_removal
        self.jinja['sequana_version'] = version

        # Common information to be filled (possibly)
        #self.data['command'] = "unset"
        self.jinja['dependencies'] =  self.get_table_dependencies('sequana').to_html()

        # the menu has a back button that may not always be the index.html
        self.jinja["main_link"] = output_filename
        self.input_filename = "undefined"

        self.jinja['online_link'] = "http://sequana.readthedocs.io/en/latest/pipelines.html"


        # Another set of data for the HTML is the galleria them
        if init_report:
            import shutil
            target = directory + "/galleria/themes"
            try:
                shutil.copytree(sequana_path + "/sequana/resources/js/galleria/themes",
                    target)
            except:
                pass
        else:
            try:
                os.makedirs(directory)
            except:
                pass

    def parse(self):
        """populate the :attr:`data` attribute used by the JINJA templates

        """
        raise NotImplementedError

    def create_report(self, onweb=False):
        try:
            self.parse()
        except Exception as err:
            raise(err)
        super(BaseReport, self).create_report(onweb=onweb)

    def read_snakefile(self, snakefile):
        """

        :param str filename:
        """
        with open(snakefile, "r") as fin:
            self.jinja['snakefile'] = fin.read()

    def read_configfile(self, configfile):
        with open(configfile, "r") as fin:
            self.jinja['config'] = fin.read()

    def copy_html_to_report(self, pattern):
        easydev.DevTools().mkdir(self.directory)
        return self._copy_files(pattern, self.directory)

    def _copy_files(self, pattern, where):
        sources = glob.glob(pattern)
        targets = []
        for source in sources:
            # ignore the report directory itself
            if source.startswith(self.directory):
                continue
            filename = source.rsplit("/", 1)[1]
            target = where + "/" +filename
            # if files are identical, fails
            try:
                shutil.copy(source, target) 
            except:pass
            targets.append(target.replace("report/", ""))
        return targets

    def copy_images_to_report(self, pattern):
        easydev.DevTools().mkdir("report")
        easydev.DevTools().mkdir("report/images")
        return self._copy_files(pattern, "report/images")

    def onweb(self):
        from easydev import onweb
        onweb(self.directory + os.sep + self.filename)

    def htmltable(self, df, tablename, bgcolors=[]):
        from reports import HTMLTable
        class SequanaHTMLTable(HTMLTable):
            def __init__(self, df, tablename, directory, bgcolors=[]):
                super(SequanaHTMLTable, self).__init__(df, tablename)
                self.directory = directory
                _basename = self.directory + os.sep + self.name
                self.csv_filename = _basename + ".csv"
                self.json_filename = _basename + ".json"
                self.bgcolors = bgcolors
            def to_html(self, index=False):
                local_html = HTMLTable(self.df.copy())
                for color in self.bgcolors:
                    local_html.add_bgcolor(color)
                html = local_html.to_html(index=index)
                pattern = '<div>%s <p>Download the CSV <a href="%s">file</a> </p> </div>' 
                return pattern % (html, self.name + ".csv")
            def to_csv(self, index=False, header=True):
                self.df.to_csv(self.csv_filename, sep=",", index=index, header=True)
            def to_json(self, filename=None):
                pass
        html = SequanaHTMLTable(df, tablename, self.directory, bgcolors)
        return html


class SequanaReport(BaseReport):
    """The main sequana report


    """
    def __init__(self, sample_dict, sample, snakefile="Snakefile",
                 configfile="config.yaml", stats="stats.txt",
                 directory="report", output_filename="index.html",  **kargs):
        """.. rubric:: constructor

        :param dict sample_dict  dictionnary of samples
        :param str sample: sample name
        :param str snakefile: the filename of the snakefile
        :param str configfile: the name of the snakemake config file
        :param str stats: where to save the stats
        :param str directory: where to save the report
        :param str output_filename: the name of the final HTML file

        This class reads the snakefile and config file to copy them in the
        directory and the HTML documents.


        """
        super(SequanaReport, self).__init__(jinja_filename="main/index.html",
                                            directory=directory,
                                            output_filename=output_filename,
                                            **kargs)

        if snakefile: 
            self.read_snakefile(snakefile)

        if configfile:
            self.read_configfile(configfile)

        try:
            self.jinja['project'] = sample 

            html = ""
            for link in dict_list[sample]:
                filename = link.split("/")[-1]
                html += '<li><a href="%s">%s</a></li>\n' % (link, filename)
                html += "</ul>"
                self.jinja['dataset'] = html
        except:
            pass

        self.jinja['snakemake_stats'] = "snakemake_stats.png"
        self.jinja['title'] = "Sequana Report"

    def parse(self):
        pass
