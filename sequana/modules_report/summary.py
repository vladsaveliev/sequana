# coding: utf-8
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
"""Module to write summary.html have all information about the pipeline and
to visit other analysis"""
import easydev
from sequana.modules_report.base_module import SequanaBaseModule

class SequanaModule(SequanaBaseModule):
    """ Write summary HTML report of an analysis. It contains all information
    about the pipeline used, input/output files and version of software.
    """
    def __init__(self, data):
        """
        """
        super().__init__()
        self.title = "Sequana Report Summary"
        self.create_report_content
        self.create_html("summary.html")
        self.json = data
        
    def create_report_content(self):
        """ Create the report content.
        """
        self.sections = list()

        self.pipeline_inputs()
        self.pipeline_outputs()
        self.workflow()
        self.running_stats()
        self.dependencies()

    def pipeline_inputs(self):
        """ Link to inputs analysed.
        """
        pass

    def pipeline_outputs(self):
        """ Link to important output generate by the pipeline
        """
        pass

    def workflow(self):
        """ Create the interactive DAG to navigate through pages.
        """
        # move the SVG file in the images directory
        img_dir = config.output_dir + os.sep + "images"
        img = self.copy_file(self.json['rulegraph'], img_dir)
        dag_svg = self.include_svg_image(img)
        with open(config.snakefile, 'r') as fp:
            code = self.add_code_section(fp.read(), 'python')
        sf = self.create_hide_section('Sf', "Show/hide Snakemake file",code)
        sf = "\n".join(sf)
        with open(config.config, 'r') as fp:
            code = self.add_code_section(fp.read(), 'yaml')
        c = self.create_hide_section('C', "Show/hide config file", code)
        c = "\n".join(c)
        self.sections.append({
            'name': 'Workflow',
            'anchor': 'workflow',
            'content': (
                "<p>The following network shows the workflow of the pipeline. "
                "Blue boxes are clickable and redirect to dedicated reports."
                "</p>\n{0}\n"
                "<p>The analysis was performed with the following Snakemake "
                "and config file:</p>\n"
                "<u><li>{1}</li>\n<li>{2}</li></u>\n".format(dag_svg, sf, c))
        })

    def running_stats(self):
        """ Barplot that shows computing time of each rule.
        """
        png = self.png_to_embedded_png(config.stats_png)
        l, c = self.create_hide_section('Stats', 'collapse/expand', png)
        self.sections.append({
            'name': "Running Stats {0}".format(self.add_float_right(l)),
            'anchor': 'stats',
            'content': c
        })

    def dependencies(self):
        """ Table with all python dependencies and a text file with tools
        needed and its versions.
        """
        html_table = self.get_table_dependencies()
        pypi = self.create_link('Pypi', 'http://pypi.python.org')
        req = self.create_link('requirements', config.requirements)
        content = ("<p>Python dependencies (<b>{0}</b>){1}</p>"
                   "<p>Dependencies downloaded from bioconda"
                   "<b>{2}</b></p>".format(pypi, html_table, req))
        l, c = self.create_hide_section('Dep', 'collapse/expand', content)
        self.sections.append({
            'name': "Dependencies {0}".format(self.add_float_right(l)),
            'anchor': 'dependencies',
            'content': c
        })

    def get_table_dependencies(self):
        """ Return dependencies of sequana.
        """
        dep_list = easydev.get_dependencies('sequana')
        version = [dep.version for dep in dep_list]
        pypi = 'https://pypi.python.org/pypi/{0}'
        url = [self.create_link(dep.project_name,
               pypi.format(dep.project_name)) for dep in dep_list]
        version = [dep.version for dep in dep_list]
        df = pd.DataFrame({'package': url, 'version': version})
        html = self.dataframe_to_html_table(df)
        return html