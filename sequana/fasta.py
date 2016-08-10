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
import textwrap
import os

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

    def format_contigs_denovo(self, project, output_dir=".", len_min=500):
        """ Method to replace NODE with the project name and to generate two
        fasta files with contigs taller than len_min and contigs smaller than
        len_min. Contigs names must be with this syntax (default syntax of 
        spades and velvet):
            NODE_1_length_524827_cov_9.49275

        :param str project: project name for output and contigs names.
        :param int cov_min: minimal length of contigs.

        Example:
        
            from sequana import FastA

            contigs = FastA("denovo_assembly.fasta")
            contigs.format_contigs_denovo(project1)

        Results are stored in files project1.ab500.fasta (above cov_min) and
        project1.bl500.fastai (below cov_min).
        """
        # check if directory exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        basename = output_dir + os.sep + project
        with open("{}.ab{}.fasta".format(basename, len_min), "w") as ab_out:
            with open("{}.bl{}.fasta".format(basename, len_min), "w") as bl_out: 
                for contigs in self:
                    name = contigs.name.split("_")
                    new_name = ">{}_{} {}\n".format(project, name[1], 
                            "_".join(name[2:]))
                    sequence = textwrap.fill(contigs.sequence, width=80) + "\n"
                    if int(name[3]) < len_min:
                        bl_out.write(new_name + sequence)
                    else:
                        ab_out.write(new_name + sequence)
