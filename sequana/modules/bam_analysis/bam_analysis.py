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

""" Sequana module to parse JSON from sequana.bamtools.
"""

class SequanaModuleReport(BaseSequanaModuleReport):

    def __init__(self):
        super().__init__(name="Sequana bam_analysis", mod_name="bam_analysis",
                         href="", info=("is a module to control the quality "
                         "of a mapping."))

