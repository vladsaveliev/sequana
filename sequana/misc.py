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


def wget(link, output):
    """Retrieve a file from internet.

    :param str link: a valid URL
    :param str output: the output filename

    .. warning:: no sanity check of any kind for now
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


def moving_average(data, n):
        """Compute moving average

        :param n: window's size.

        """
        ret = np.cumsum(data, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ma = ret[n - 1:] / n
        return ma


class FixFailedProjects(object):
    """Search for projects not finished and create a new script

    ::

        fixer = FixFailedProjects()
        # executes the method create_new_starter() and save in multirun.sh
        fixer()

    :Motivation: When running a snakemake on a large set of fastq pairs, 
        one should use the sequana utility with the --glob option. This creates 
        as many directories as pairs of fastq files. each directory with a
        snakefile and a runme.sh script. Now that's a lot of projects so to
        enter into each directory and run the scripts, we also have a 
        general script  **multirun.sh** that contains all commands. Once ran,
        some jobs may fail. This is now difficult to know which one have failed.
        This class will figure it out and create a new multirun.sh file

    """
    def __init__(self, where=".", pattern="17 of 17"):
        
        self.pattern = pattern
        self.where = where
        self._template = """
        cd %(directory)s
        sh runme.sh &
        sleep 0.5
        echo Starting %(directory)s
        cd ..
        """

    def __call__(self):
        res = self.create_new_starter()
        # figure out the new file


    def create_new_starter(self):
        count = 0
        txt = ""
        for directory in glob.glob("*"):
            if os.path.isdir(directory) is False:
                continue

            # else:
            filename = directory + "/run.err"; 
            if os.path.exists(filename) is False:
                txt += self._template % {"directory":directory}
                count += 1
            else:
                data = open(filename, "r").read()
                if self.pattern not in data:
                    txt += self._template % {"directory":directory}
                    count += 1
        return txt

