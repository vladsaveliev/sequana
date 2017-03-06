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
""" Utils to create a Jquery DataTables for your HTML file. First class
contains the javascript function to create DataTable. The second set table
which use javascript created in the first class.
"""
from collections import OrderedDict

from sequana import logger


class DataTableFunction(object):
    """ Class that contains Jquery DataTables function and options.

    Example:

    ::

        df = pandas.read_csv('data.csv')
        datatable_js = DataTableFunction(df, 'data')
        datatable_js.datatable_options = {'pageLength': 15,
                                          'dom': 'Bfrtip',
                                          'buttons': ['copy', 'csv']}
        js = datatable_js.create_javascript_function()
        html_datatables = [DataTable(df, "data_{0}".format(i), datatable_js)
                           for i, df in enumerate(df_list)]

    All options of datatable:
        https://datatables.net/reference/option/
    """
    def __init__(self, df, html_id):
        """.. rubric:: contructor

        :param df: data frame.
        :param str html_id: the ID used in the HTML file.
        """
        self._html_id = html_id
        self._datatable_options = dict()
        self._datatable_columns = self._set_datatable_columns(df)

    @property
    def html_id(self):
        """ Get the html_id. Html_id cannot be set by the user after the
        initiation of the class.
        """
        return self._html_id

    @property
    def datatable_options(self):
        """ Get, set or delete the DataTable options. Setter take a dict as
        parameters with desired options and update the dictionary.

        Example:

        ..
            datatable = DataTableFunction("tab")
            datatable.datatable_options = {'dom': 'Bfrtip',
                                           'buttons': ['copy', 'csv']}

        source: https://datatables.net/reference/option/
        """
        return self._datatable_options

    @datatable_options.setter
    def datatable_options(self, d):
        self._datatable_options.update(d)

    @datatable_options.deleter
    def datatable_options(self):
        self._datatable_options = dict()

    @property
    def datatable_columns(self):
        """ Get datatable_columns dictionary. It is automatically set from the
        df you want plot.
        """
        return self._datatable_columns

    def _set_datatable_columns(self, df):
        """ Fill :attr:`DataTableFunction.datatable_columns` with header of
        :param:`DataTableFunction.df`.
        """
        column_dict = OrderedDict((name, dict()) for name in df.columns)
        return column_dict

    def create_javascript_function(self):
        """ Return javascript to create the DataTable.
        """
        js_function = """
<script type="text/javascript">
    function parseCsv_{0}(csv, id) {{
        Papa.parse(csv, {{
            comments: '#',
            delimiter: ',',
            header: true,
            dynamicTyping: true,
            error: function(reason) {{
                console.log(reason);
            }},
            complete: function(results) {{
                {1}
            }}
        }});
    }};
</script>
"""
        return js_function.format(self.html_id,
                                  self._create_datatable_option())

    def _create_datatable_option(self):
        """ Return DataTable options.
        """
        self.datatable_options['columns'] = self._create_columns_option()
        js = self._dict_to_string(self.datatable_options)
        js = "$(id).DataTable({{{0},data: results.data}});".format(js)
        return js

    def _create_columns_option(self):
        """ Return string well formated with all columns options.
        """
        js = [self._coloption_2_str(key, value) for key, value in
              self.datatable_columns.items()]
        return '[{0}]'.format(',\n'.join(js))

    def _coloption_2_str(self, name, options):
        s = "data:'{0}'".format(name)
        if options:
            s = "{0},\n{1}".format(s, self._dict_to_string(options))
        return '{{{0}}}'.format(s)

    def _dict_to_string(self, d):
        """ Convert dict to string for CanvasJS.

        Example:

        ::
            dico = {'key1': value1, 'key2': value2, 'key3': value3}
            print(CanvasJS._dict_to_string(dico))

            "key1:value1,key2:value2,key3:value3"
        """
        s = ['{0}:{1}'.format(key, self._check_type(value)) for key, value in
             d.items()]
        return ',\n'.join(s)

    def _check_type(self, value):
        """ Check type of value to fill javascript sections. String must be
        surrounded by quote and not boolean or integer.

        Javascript variable must not be surrounded by quote. Custom variables
        start with 'data_'.
        """
        try:
            if not value.startswith(('true', 'false', 'function', '{', '[')):
                return "'{0}'".format(value)
        except AttributeError:
            return value
        return value

    def set_links_to_column(self, link_col, target_col):
        """ Hidden a column with url and connect its with some columns.

        :param: str link_col: column with your URL.
        :param: str target_col: column to connect.
        """
        # hidden the link column
        try:
            self.datatable_columns[link_col]['visible'] = 'false'
        except KeyError:
            logger.warning("KeyError: Column name '{0}' does not exist."
                           .format(target_col))
            pass
        # function to add link
        fct = """function(data, type, row, meta){{
            return '<a href="'+row.{0}+'" target="_blank">'+data+'</a>';
        }}
        """.format(link_col)
        try:
            self.datatable_columns[target_col]['render'] = fct
        except KeyError:
            logger.warning("KeyError: Column name '{0}' does not exist."
                           .format(target_col))
            pass


class DataTable(object):
    """ Class that contains html table which used a javascript function. You
    must add in your HTML file the js function
    (:meth:`DataTable.create_javascript_function`) and the HTML code
    (:meth:`DataTable.create_datatable`).

    Example:

    ::

        df = pandas.read_csv('data.csv')
        datatable = DataTable(df, 'data')
        datatable.datatable.datatable_options = {'pageLength': 15,
                                                 'dom': 'Bfrtip',
                                                 'buttons': ['copy', 'csv']}
        js = datatable.create_javascript_function()
        html = datatable.create_datatable()
        # Second CSV file with same format
        df2 = pandas.read_csv('data2.csv')
        datatable2 = DataTable(df2, 'data2', datatable.datatable)
        html2 = datatable.create_datatable()
    """
    def __init__(self, df, html_id, datatable=None):
        """.. rubric:: contructor

        :param df: data frame.
        :param str html_id: the ID used in the HTML file.
        :param DataTableFunction datatable: javascript function to create the
            Jquery Datatables. If None, a :class:`DataTableFunction` is
            generated from the df.
        """
        self._df = df
        self._html_id = html_id
        if datatable:
            self.datatable = datatable
        else:
            self.datatable = DataTableFunction(df, html_id)

    def __len__(self):
        return len(self.df)

    @property
    def df(self):
        return self._df

    @property
    def html_id(self):
        return self._html_id

    def create_datatable(self, style="width:100%", **kwargs):
        """ Return string well formated to include in a HTML page.

        :param str style: css option of your table.
        :param **dict kwargs: parameters of :meth:`pandas.DataFrame.to_csv`.
        """
        html = """
<script type="text/javascript">
    $(document).ready(function() {{
        var {0} = document.getElementById('csv_{0}').innerText;
        parseCsv_{1}({0}, '#table_{0}');
        {0} = null;
    }});
</script>
        """.format(self.html_id, self.datatable.html_id)
        html += self._create_hidden_csv(**kwargs)
        html += self._create_html_table(style)
        return html

    def _create_hidden_csv(self, **kwargs):
        """ Return the HTML code and the CSV code for your hidden CSV section.

        :param **dict kwargs: parameters of :meth:`pandas.DataFrame.to_csv`.
        """
        csv = self._df.to_csv(**kwargs)
        html = '<pre id="csv_{0}">{1}</pre>'.format(self.html_id, csv.strip())
        css = '<style>#csv_{0}{{display:none}}</style>'.format(self.html_id)
        return '{0}\n{1}\n'.format(css, html)

    def _create_html_table(self, style):
        """ Just for set some option and header.

        :param str style: css option of your table.
        """
        # set table id
        if style:
            style = 'style="{0}"'.format(style)
        html_table = '<table id="table_{0}" class="display datatable" {1}>'\
            .format(self.html_id, style)
        # create table's header
        th = '<th>{0}</th>'
        header = [th.format(name) for name in self.df]
        header = '<thead>{0}</thead>'.format("\n".join(header))
        html_table = """
    {0}
        {1}
    </table>
        """.format(html_table, header)
        return html_table

    def create_javascript_function(self):
        """ Generate the javascript function to create the DataTable in a HTML
        page.
        """
        return self.datatable.create_javascript_function()
