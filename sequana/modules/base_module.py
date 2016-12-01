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

""" Sequana module base class contains JSON parser and some generics
functions.
"""

from sequana.utils import report


class BaseSequanaModuleReport(object):

    def __init__(self, name='base', mod_name='base', target='', href='',
                 info='', extra=''):
        self.name = name
        self.mod_name = mod_name
        if not target:
            target = self.name
        self.intro = ('<p><a href="{0}" target="_blank">{1}</a> '
                      '{2}</p>{3}').format(href, target, info, extra)
        
    def get_json_files(self, module):
        """ Search the analysis directory for JSON files produce by Sequana.
        JSON file must containt the module field to know its origins.

        :param string module: sequana module name
        """
        for json_fl in report.json_list:
            if json_fl["module"] == self.mod_name:
                yield json_fl
