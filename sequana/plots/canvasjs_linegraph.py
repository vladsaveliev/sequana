# coding: utf-8
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
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
"""Sequana class to plot a CanvasJS linegraph from an embedded csv file.
"""
from sequana import logger
from sequana.plots.canvasjs_base import CanvasJS


"""
- Need an embedded csv file with a specified id. The id is necessary to
recover data. (csvdata for example).
- Create the function to parse this specific csv file:
    - Array needed to create (cov, mapq0, gc content)
    - Function to generate canvasJS depend of previous array (add data - change
        title)

Maybe do a class with necessary code and do method to create line in plot
according to data present in csv file.

Entry:
    - Need csv file in a pd.DataFrame:
        We consider the first column is the x and other columns are y to plot.
        The plot can have two Y axis.
    - id name of the hidden csv in html page
        Id is the important to know where are datas.

Required parameter for some function:
    - which column you want plot.
        To plot a column you need to set an Y-axis setup.

Return:
    - Javascript function
    - Javascript script to execute function when document is ready
    - Hidden section with data in csv format
"""


class CanvasJSLineGraph(CanvasJS):
    """ Class to create a CanvasJS linegraphe for an HTML page. It creates a
    hidden pre section with your CSV. It is necessary embedded because browsers
    forbid the reading of data stored locally. Your html page need CanvasJS and
    PapaParse.
    """

    def __init__(self, csv, html_id, x_column, y_columns):
        """.. rubric:: constructor

        :param str csv: data as CSV format.
        :param str html_id: the ID used in your html. All function
            will have this tag.
        :param str x_column: column used as axis X for plot.
        :param list y_columns: colums used as axis Y for plot.
        """
        super().__init__(html_id)
        self.csv = csv
        self.html_id = html_id
        self.x_column = x_column
        self.y_columns = y_columns

    def create_hidden_csv(self):
        """ Return the HTML code and the CSV code for your hidden CSV section.
        """
        html = '<pre id="{0}">{1}</pre>'.format(self.html_id, self.csv)
        css = '#{0}{{display:none}}'.format(self.html_id)
        return html, css

    def create_js_csv_parser(self):
        """ Create the efficient javascript csv parser with PapaParse.

        :param str x_column: name of your data for X axis.
        :param list y_columns: names of your data you want for Y axis.
        """
        # Create variable name
        variable = ["data_{0}".format(name) for name in self.y_columns]
        # Create variable for javascript array
        init_var = ["var {0} = [];".format(name) for name in variable]
        init_var = "\n".join(init_var)
        # Fill arrays
        fill_array = ["{0}.push({{x: curow.{1}, y: curow.{2}}});".format(v,
                      self.x_column, n) for n, v in
                      zip(self.y_columns, variable)]
        fill_array = "\n".join(fill_array)
        self.variables = ", ".join(variable)
        function = """
    function processData_{0}(csv) {{
        {1}
        var curow;
        Papa.parse(csv, {{
            comments: '#',
            delimiter: ',',
            header: true,
            dynamicTyping: true,
            error: function(reason) {{
                console.log(reason);
            }},
            step: function(row){{
                curow = row.data[0];
                {2}
            }},
            complete: function(results) {{
                drawChart_{0}({3});
            }}
        }});
    }};
        """.format(self.html_id, init_var, fill_array, self.variables)
        logger.info("There are {0} variables created.".format(len(variable)))
        self.data_section = [{'dataPoints': var} for var in variable]
        return function

    def set_x_axis(self, axis_attr=dict()):
        """ Method to configure X axis of the line graph.

        :param dict axis_attr: dictionary with canvasjs axisX
            Attributes.

        axisX: http://canvasjs.com/docs/charts/chart-options/axisx/

        Example:

        ::

            line_graph = CanvaJSLineGraph(csv, csvdata)
            axisX_section = line_graph.set_x_axis({
                'title': 'Title of X data',
                'titleFontSize': 16,
                'labelFontSize': 12,
                'lineColor': '#5BC0DE',
                'titleFontColor': '#5BC0DE',
                'labelFontColor': '#5BC0DE'})
        """
        self._set_axis("axisX", axis_attr)

    def set_y_axis(self, axis_attr=dict()):
        """ Method to configure first axis Y of the line graph.

        :param dict axis_attr: dictionary with canvasjs axisY
            Attributes.

        axisY: http://canvasjs.com/docs/charts/chart-options/axisy/

        Example:

        ::

            line_graph = CanvaJSLineGraph(csv, csvdata)
            line_graph.set_y_axis({'title': 'Title of Y data',
                                   'titleFontSize': 16,
                                   'labelFontSize': 12,
                                   'lineColor': '#5BC0DE',
                                   'titleFontColor': '#5BC0DE',
                                   'labelFontColor': '#5BC0DE'})
        """
        self._set_axis('axisY', axis_attr)

    def set_y2_axis(self, axis_attr=dict()):
        """ Method to configure second axis Y of the line graph.

        :param dict axis_attr: dictionary with canvasjs axisY
            Attributes.

        axisY: http://canvasjs.com/docs/charts/chart-options/axisy/

        Example:

        ::

            line_graph = CanvaJSLineGraph(csv, csvdata)
            axisX_section = line_graph.set_x_axis({
                'title': 'Title of Y data',
                'titleFontSize': 16,
                'labelFontSize': 12,
                'lineColor': '#5BC0DE',
                'titleFontColor': '#5BC0DE',
                'labelFontColor': '#5BC0DE'})
        """
        self._set_axis('axisY2', axis_attr)

    def create_canvas_js_object(self):
        """ Method to convert all section as javascript function.

        Return command line to launch generation of plot, js function to create
        CanvasJS object and the html div that contains CanvasJS plot.
        """
        # js function to create the chart container
        js = """
    function drawChart_{0}({1}) {{
        console.log("Start drawChart");
        var chart = new CanvasJS.Chart("chartContainer_{0}",{2}
        )
        chart.render()
        }};
        """
        canvas_attr = super().create_canvas_js_object()
        canvas_js = js.format(self.html_id, self.variables, canvas_attr)
        # js command to run the CanvasJS
        js_command = """
        var csv_{0} = document.getElementById('{0}').innerText;
        processData_{0}(csv_{0});
        """.format(self.html_id)
        canvas_html = ('<div id="chartContainer_{0}" style="height: '
                       '450px; width: 100%;"></div>'.format(self.html_id))
        return js_command, canvas_js, canvas_html
