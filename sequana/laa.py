# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################




from sequana import BAM
import glob
import pandas as pd
import pylab
from sequana import logger
logger.name = __name__


class LAA():
    def __init__(self, where="bc*"):
        self.filenames = glob.glob(where + "/" + "amplicon_*summary.csv")
        self.data = [pd.read_csv(this) for this in self.filenames]
        self.numbers = [len(x) for x in self.data]
        self.df = pd.DataFrame(
            [(x[0:4],y) for x, y in zip(self.filenames, self.numbers)])

    def hist_amplicon(self, fontsize=12):
        pylab.hist(self.numbers, bins=max(data), ec="k", align="left")
        pylab.ylabel("#", fontsize=fontsize)
        pylab.ylabel("Number of amplicons per barcode", fontsize=fontsize)

    def __str__(self):
        return self.df.sort_values(1).to_string()

    def plot_max_length_amplicon_per_barcode(self):
        data = [max(x.SequenceLength) if len(x.SequenceLength) else 0  for x in self.data]
        pylab.plot([int(x[2:4]) for x in self.filenames], data, "o")
        pylab.xlabel("barcode name (filename)")
        pylab.ylabel("Max length amongst amplicons (0 if no amplicon)")
        #savefig("max_amplicon_length_per_barcode.png", dpi=200)


class LAA_Assembly():
    """

    Input is a SAM/BAM from the mapping of amplicon onto a known reference.
    Based on the position, we can construct the new reference.

    """
    def __init__(self, filename):
        self.bam = BAM(filename)


    def build_reference(self):
        self.bam.reset()
        # scan BAM file assuming it is small
        aa = [a for a in self.bam]

        # retrieve data of interest
        data = [(a.pos, {
                    "name":a.query_name,
                    "sequence": a.query_sequence,
                    "cigar": a.cigarstring,
                    "position": a.pos,
                    "qstart": a.qstart,
                    "qend": a.qend}) for a in aa]

        # sort by starting position
        data.sort(key=lambda x: x[0])

        for i, read in enumerate(data):
            read = read[1]
            if i == 0:
                sequence = read["sequence"]     # 2 is query_sequence
            else:
                pr = data[i-1][1]   # previous read
                L = len(pr["sequence"])
                end_position_pr = pr['position'] - pr['qstart'] + L 

                # overlap between previous read and this one
                overlap = end_position_pr - (read['position'] - read['qstart']) +0
                print(overlap)
                print(pr['position'], pr['qstart'], L, end_position_pr)
                print(read['position'], read['qstart'])
                sequence = sequence + read["sequence"][overlap+1:]

        # argmax([sum(a==b for a,b in zip(X[-i:] , Y[:i]))/float(i+1) for i in range(1000)])
        return sequence

    def save_fasta(self, filename, sequence=None):
        if sequence is None:
            sequence = self.build_reference()


        with open(filename, "w") as fout:
            fout.write(">test\n{}".format(sequence))


