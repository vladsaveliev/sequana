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
"""misc utilities"""
import os
import glob
import numpy as np


def textwrap(text, width=80, indent=0):
    if indent == 0:
        indent = ""
    else:
        indent = " " * indent
    data = [indent + text[i*width:(i+1)*width:] for i in range(len(text)//width + 1)]
    return "\n".join(data)


def wget(link, output):
    """Retrieve a file from internet.

    :param str link: a valid URL
    :param str output: the output filename

    .. warning:: no sanity check of any kind for now
    .. todo:: move to easydev
    """
    try:
        from urllib import urlretrieve
    except:
        from urllib.request import urlretrieve
    urlretrieve(link, filename=output)


def findpos(seq, chr):
    N = len(chr)
    for i, dummy in enumerate(seq):
        if seq[i:i+N] == chr:
            yield i

