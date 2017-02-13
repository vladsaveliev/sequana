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




def galleria(filenames):

    html = """ <div class="galleria">"""
    for filename in filenames:
        html += """<img src="%s" 
                        data-title="%s" 
                        data-description="%s">\n""" % (filename, filename, filename)
    html += "</div>"
    html += """
 <script>
     Galleria.loadTheme('galleria/themes/classic/galleria.classic.min.js');
     Galleria.run('.galleria');
 </script>
    """

    return html

