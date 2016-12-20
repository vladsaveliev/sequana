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
"""Generic module is the parent module of all other module"""
import jinja2
import io
import base64

from reports import HTMLTable


class GenericModule(object):
    """ Generic Module to write HTML reports.
    """
    def __init__(self,
                 template,
                 output_directory="report/"):
        """
        """
        env = jinja2.Environment(loader=jinja2.FileSystemLoader(template))
        self.j_template = env.get_template("base.html")

    def dataframe_to_html_table(self, dataframe, index):
        """
        """
        html = HTMLTable(dataframe)
        return html.to_html(index=index)

    def create_embed_png(self, plot_function, kwargs=dict(), style=None):
        """ Take as a plot function as input and create a html embed png image.
        """
        buf = io.BytesIO()
        plot_function(filename=buf, **kwargs)
        if style:
            html = '<img style="{0}" '.format(style)
        else:
            html = '<img '
        html += 'src="data:image/png;base64,{0}"/>'.format(
            base64.b64encode(buf.getvalue()).decode("utf-8"))
        buf.close()
        return html
