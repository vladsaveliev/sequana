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
""" Sequana report config contains configuration informations to create HTML
report with Jinja2.
"""
import sys
import os
import pkg_resources
import glob

import easydev


# Set logger
from sequana import logger

# Get sequana informations
version = pkg_resources.get_distribution('sequana').version
script_path = os.path.dirname(os.path.realpath(__file__))

# Default options
output_dir = os.path.realpath(os.getcwd())
template = "standard"

def get_entry_name(entry_point):
    return str(entry_point).split('=')[0].strip()

# Templates available
template_dict = {get_entry_name(entry_point): entry_point for entry_point
                 in pkg_resources.iter_entry_points('sequana.report_template')}

# Modules available
module_dict = {get_entry_name(entry_point): entry_point for entry_point
                 in pkg_resources.iter_entry_points('sequana.module')}

# Check if templates are found
if len(template_dict) == 0:
    print("Error - Sequana report templates are not found", file=sys.stderr)
    sys.exit(1)

# Find css/js file necessary for the report
sequana_path = easydev.get_package_location('sequana')
css_path = os.sep.join([sequana_path, "sequana", "resources", "css"])
js_path = os.sep.join([sequana_path, "sequana", "resources", "js"])
css_list = glob.glob(css_path + os.sep + "*css")
js_list = glob.glob(js_path + os.sep + "*js")

# Sections list for summary.html
summary_sections = list()
