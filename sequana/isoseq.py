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
import glob
import os

from sequana.lazy import  pandas as pd
from sequana import FastQ, FastA
from sequana.lazy import pylab


__all__ = ["IsoSeq"]


class IsoSeq(object):
    """


    ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/

    """
    def __init__(self, directory=".", prefix="job-*"):
        self.prefix = prefix
        self.directory = directory

        # low quality isoforms
        self.lq_isoforms = self.get_file("lq_isoforms.fastq")
        if self.lq_isoforms:
            self.lq_sequence = FastQ(self.lq_isoforms)

        # high quality isoforms
        self.hq_isoforms = self.get_file("hq_isoforms.fastq")
        if self.hq_isoforms:
            self.hq_sequence = FastQ(self.hq_isoforms)

        # General info
        self.csv = self.get_file("-file.csv")
        if self.csv:
            self.data = pd.read_csv(self.csv)

        # CCS fasta sequence
        #self.ccs = self.get_file("-ccs.tar.gz")
        self.ccs = self.get_file("ccs.fasta", noprefix=True)
        if self.ccs:
            self.ccs = FastA(self.ccs)

    def get_file(self, tag, noprefix=False):
        if noprefix:
            filenames = glob.glob(self.directory + os.sep + tag)
        else:
            filenames = glob.glob(self.directory + os.sep + self.prefix + tag)
        if len(filenames) == 1:
            return filenames[0]
        elif len(filenames) > 1:
            print("Found several files ending in %s" % tag)
        else:
            print("No files matching %s" % tag)
        return None

    def stats(self):
        results = {}
        if self.data is not None:
            print("Reading strand")
            results['strand'] = {
                "+": sum(self.data.strand == "+"),
                "-": sum(self.data.strand == "-"),
                "?": sum(self.data.strand.isnull())
            }

            results['classification'] = {
                "total_ccs_reads" : len(self.data),
                "five_prime_reads" : sum(self.data.fiveseen),
                "three_prime_reads" : sum(self.data.threeseen),
                "polyA_reads" : sum(self.data.polyAseen),
            }

        if self.lq_isoforms:
            print("Reading LQ isoform")
            results['lq_isoform'] = self.lq_sequence.stats() # number of 

        if self.hq_isoforms:
            print("Reading HQ isoform")
            results['hq_isoform'] = self.hq_sequence.stats() # number of polished HQ isoform

        if self.ccs:
            seq = [ read.sequence for read in self.ccs]
            results["CCS"] = {
                "mean_length" : pylab.mean(seq),
                "number_ccs_bases" : sum(seq),
                "number_ccs_reads" : len(seq)
            }

        self.idents_v = []
        self.full_v = []
        self.non_full_v = []
        self.isoform_lengths = []
        for read in self.lq_sequence:
            ident, full, non_full, length = read['identifier'].decode().split(";")
            self.idents_v.append(ident)
            self.full_v.append(int(full.split("=")[1]))
            self.non_full_v.append(int(non_full.split("=")[1]))
            self.isoform_lengths.append(int(length.split("=")[1]))

        return results

    def hist_read_length_consensus_isoform(self, mode="all", bins=80, rwidth=0.8,
        align="left", fontsize=16, edgecolor="k", **kwargs):
        """

        mode can be all, lq, hq
        """

        pylab.clf()

        L1 = [len(read['sequence']) for read in self.lq_sequence]
        L2 = [len(read['sequence']) for read in self.hq_sequence]
        if mode == "all":
            L = L1 + L2
        elif mode == "lq":
            L = L1
        else:
            L = L2
 
        Y, X, _ = pylab.hist(L, bins=bins, rwidth=rwidth, align=align,
                    ec=edgecolor, **kwargs)
        pylab.gca().set_ylim(bottom=0)
        pylab.gca().set_xlim(left=0)
        pylab.xlabel("Read length", fontsize=fontsize)
        pylab.ylabel("Number of reads", fontsize=fontsize)
        pylab.grid()

        ax_twin = pylab.gca().twinx()

        shift = (X[1] - X[0]) / 2

        ax_twin.plot(X[1:]- shift, len(L)-pylab.cumsum(Y), "k")
        ax_twin.set_ylim(bottom=0)
        pylab.ylabel("CDF", fontsize=fontsize)

        pylab.title("Read length of Consensus isoforms reads")

    def hist_average_quality(self, fontsize=16):

        hq_qv = [mean([phred.ascii_to_quality(X) for X in read['quality'].decode()]) 
                for read in iso.hq_sequence]
        lq_qv = [mean([phred.ascii_to_quality(X) for X in read['quality'].decode()]) 
            for read in iso.lq_sequence]

        Y1, X = numpy.histogram(hq_qv, bins=range(0,94))
        Y2, X = numpy.histogram(lq_qv, bins=range(0,94))
        pylab.bar(X[1:], Y1, width=1, label="HQ")
        pylab.bar(X[1:], Y2, bottom=Y1, width=1, label="LQ")
        pylab.xlabel("Isoform average QV")
        pylab.ylabel("# Isoform")
        pylab.legend(fontsize=fontsize)


#other files:
"""
job-288-57b26283-8c10-71c9-d522-8397532135ff-hq_isoforms.contigset.xml
job-288-5ee91702-dfa2-2b43-8b5b-996b597ed650-file.consensusreadset.xml
job-288-c6d3ccef-afdf-aa13-f55f-30b8b4d8a375-lq_isoforms.contigset.xml
"""

