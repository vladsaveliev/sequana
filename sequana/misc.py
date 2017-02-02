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
""".. rubric:: misc utilities"""
import os
import glob
import numpy as np
import platform


def textwrap(text, width=80, indent=0):
    """Wrap a string with 80 characters

    :param text: input text
    :param width: (defaults to 80 characters)
    :param indent: possible indentation (0 by default)

    """
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
    """Find position(s) of a substring into a longer string.

    Note that this function is a generator::

        >>> list(findpos("AACCGGAAGGTT", "GG"))
        [4,8]

    """
    N = len(chr)
    for i, dummy in enumerate(seq):
        if seq[i:i+N] == chr:
            yield i


def on_cluster(pattern=["tars-"]):
    """Used to check if we are on a cluster

    "tars-" is the name of a cluster's hostname.
    Change or append the argument **pattern** with your cluster's hostname

    :param str pattern: a list of names (strings) or a string

    """
    if isinstance(pattern, str):
        pattern = [pattern]

    for this in pattern:
        if platform.uname().node.startswith(this):
            return True
        else:
            return False
