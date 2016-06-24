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
"""Utilities to manipulate FASTQ and Reads

"""
from pysam import FastxFile


__all__ = ["FastA"]

# cannot inherit from FastxFile (no object in the API ?)
class FastA(object):
    """Class to handle FastA files


    """
    def __init__(self, filename, verbose=False):
        self._fasta = FastxFile(filename)
        self._N = len([x for x in FastxFile(filename)])

    def __iter__(self):
        return self

    def __next__(self): # python 3
        return self.next()

    def next(self): # python 2
        # reads 4 lines
        try:
            d = next(self._fasta)
            return d
        except KeyboardInterrupt:
            # This should allow developers to break a loop that takes too long
            # through the reads to run forever
            self._fasta.close()
            self._fasta = FastxFile(self._fasta.filename)
        except:
            self._fasta.close()
            self._fasta = FastxFile(self._fasta.filename)
            raise StopIteration
        return d

    def __len__(self):
        return self._N

    def _get_names(self):
        return [this.name for this in self]
    names = property(_get_names)

    def _get_sequences(self):
        return [this.sequence for this in self]
    sequences = property(_get_sequences)

    def _get_comment(self):
        return [this.comment for this in self]
    comment = property(_get_comment)





