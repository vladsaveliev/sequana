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
import easydev
sepjoin = os.sep.join

from reports import Report
from sequana import version

import glob


class BaseReport(Report):
    """A Parent child for all reports created in Sequana


    """
    def __init__(self, jinja_filename, directory="report",
                output_filename="test.html", **kargs):
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
            **kargs)

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

        # Another set of data for the HTML is the galleria them
        import shutil
        target = directory + "/galleria/themes"
        try:
            shutil.copytree(sequana_path + "/sequana/resources/js/galleria/themes",
                target)
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



class SequanaReport(BaseReport):
    """The main sequana report


    """
    def __init__(self, snakefile="Snakefile", configfile="config.yaml",
                 stats="stats.txt", directory="report",
                output_filename="index.html",  **kargs):
        """.. rubric:: constructor


        :param str snakefile: the filename of the snakefile
        :param str configfile: the name of the snakemake config file
        :param str stats: where to save the stats
        :param str directory: where to save the report
        :param str output_filename: the name of the final HTML file

        This class reads the snakefile and config file to copy them in the
        directory and the HTML documents.


        """

        super(SequanaReport, self).__init__(
            jinja_filename="main/index.html",
            directory=directory,
            output_filename=output_filename,
            **kargs)

        self.read_snakefile(snakefile)
        self.read_configfile(configfile)

        try:
            from sequana.snaketools import SequanaConfig
            config = SequanaConfig(configfile)
            self.jinja['project'] = config.PROJECT

            html = ""
            for link, filename in zip(config.DATASET, config.BASENAME):
                html += '<li><a href="../%s">%s</a></li>\n' % (link, filename)
                html += "</ul>"
                self.jinja['dataset'] = html
        except:
            pass

        self.jinja['online_link'] = "http://sequana.readthedocs.io/en/latest/pipelines.html"
        self.jinja['snakemake_stats'] = "snakemake_stats.png"
        self.jinja['title'] = "Sequana Report"

    def parse(self):
        pass











