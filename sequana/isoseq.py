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
from collections import Counter


from sequana.lazy import  pandas as pd
from sequana import FastQ, FastA, BAMPacbio
from sequana.lazy import pylab
from sequana import logger, phred

import numpy as np

logger.name = __name__


__all__ = ["IsoSeqQC", "IsoSeqBAM"]


class IsoSeqQC(object):
    """


    Use get_isoseq_files on smrtlink to get the proper files

    iso = IsoSeqQC()
    iso.hist_read_length_consensus_isoform() # histo CCS 
    iso.stats() # "CCS" key is equivalent to summary metrics in CCS report


    todo: get CCS passes histogram . Where to get the info of passes ? 

    """
    def __init__(self, directory=".", prefix=""):
        self.prefix = prefix
        self.directory = directory
        self.sample_name = "undefined"

        # low quality isoforms
        filename = "all.polished_lq.fastq"
        self.lq_isoforms = self.get_file(filename)
        if self.lq_isoforms:
            logger.info("Reading {}".format(filename))
            self.lq_sequence = FastQ(self.lq_isoforms)

        # high quality isoforms
        filename = "all.polished_hq.fastq"
        self.hq_isoforms = self.get_file(filename)
        if self.hq_isoforms:
            logger.info("Reading {}".format(filename))
            self.hq_sequence = FastQ(self.hq_isoforms)

        # General info
        filename = "file.csv"
        self.csv = self.get_file(filename)
        if self.csv:
            logger.info("Reading {}".format(filename))
            self.data = pd.read_csv(self.csv)

        # CCS fasta sequence
        #self.ccs = self.get_file("-ccs.tar.gz")
        filename = "ccs.fasta"
        self.ccs = self.get_file(filename, noprefix=True)
        if self.ccs:
            logger.info("Reading {}".format(filename))
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
            logger.info("Reading strand")
            results['strand'] = {
                "+": sum(self.data.strand == "+"),
                "-": sum(self.data.strand == "-"),
                "?": sum(self.data.strand.isnull())
            }

            results['classification'] = {
                "total_ccs_reads" : len(self.data),
                "five_prime_reads" : int(self.data.fiveseen.sum()),
                "three_prime_reads" : int(self.data.threeseen.sum()),
                "chimera" : int(self.data.chimera.sum()),
                "polyA_reads" : int(self.data.polyAseen.sum()),
            }

        if self.lq_isoforms:
            logger.info("Reading LQ isoforms")
            results['lq_isoform'] = self.lq_sequence.stats() # number of 

        if self.hq_isoforms:
            logger.info("Reading HQ isoforms")
            results['hq_isoform'] = self.hq_sequence.stats() # number of polished HQ isoform

        if self.ccs:
            seq = [ len(read.sequence) for read in self.ccs]
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


    def to_summary(self, filename="sequana_summary_isoseq.json", data=None):
        """Save statistics into a JSON file

        :param filename:
        :param data: dictionary to save. If not provided, use :meth:`stats`

        """
        from sequana.summary import Summary
        if data is None:
            data = self.stats()
        Summary("isoseq",self.sample_name, data=data).to_json(filename)


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

        ax_twin.plot(X[0:-1]- shift, len(L) - pylab.cumsum(Y), "k")
        ax_twin.set_ylim(bottom=0)
        pylab.ylabel("CDF", fontsize=fontsize)

        pylab.title("Read length of Consensus isoforms reads")

    def hist_average_quality(self, fontsize=16, bins=None):
        """

        bins is from 0 to 94 
        """

        hq_qv = [pylab.mean([ord(X)-33 for X in read['quality'].decode()]) 
                for read in self.hq_sequence]
        lq_qv = [pylab.mean([ord(X) -33 for X in read['quality'].decode()]) 
            for read in self.lq_sequence]

        if bins is None:
            bins = range(0,94)
        Y1, X = np.histogram(hq_qv, bins=bins)
        Y2, X = np.histogram(lq_qv, bins=bins)
        pylab.bar(X[1:], Y1, width=1, label="HQ")
        pylab.bar(X[1:], Y2, bottom=Y1, width=1, label="LQ")
        pylab.xlim([0.5, 93.5])

        pylab.xlabel("Isoform average QV")
        pylab.ylabel("# Isoform")
        pylab.legend(fontsize=fontsize)

        ax = pylab.twinx()
        N = np.sum(Y1+Y2)
        ax.plot(X, [N] + list(N-np.cumsum(Y1+Y2)), "k")


class IsoSeqBAM(object):
    """Here, we load a BAM file generated with minimap using 
    as input the BAM file  created with the mapping og HQ isoforms
    on a reference.

    df contains a dataframe for each read found in the SAM (and hq_isoform)
    we populate the GC content, the mapping flag, the reference name (-1 means
    no mapping i.e flag ==4). flag of 4 means unmapped and there is no 
    ambiguity about it.

    In the data file example, other falgs are 0, 16 (SEQ being reverse
    complement<F12>) , 2048 (supplementary     segment).

    Example of minimap2 command::

        minimap2 -t 4  -ax splice -uf --secondary=no  SIRV-E0.fa 
            hq_isoforms.fasta > hq_isoforms.fasta.sam 2> hq_isoforms.fasta.sam.log



    """
    def __init__(self, filename):
        self.filename = filename
        self.bam = BAMPacbio(self.filename)

    @property
    def df(self):
        self.bam.reset()
        rnames = [self.bam.data.get_reference_name(a.rname) if a.rname!=-1 
                  else -1 for a in self.bam.data]

        df = self.bam.df.copy()
        df['reference_name'] = rnames

        self.bam.reset()
        df['flags'] = [a.flag for a in self.bam.data]

        self.bam.reset()
        df['mapq'] = [a.mapq for a in self.bam.data]

        self.bam.reset()
        df['cigar'] = [a.cigarstring for a in self.bam.data]

        # Drop SNR that are not populated in the mapped BAM file.
        df.drop(['snr_A', 'snr_C', 'snr_G', 'snr_T'], axis=1, inplace=True)

        return df


    def hist_isoform_length_mapped_vs_unmapped(self, bins=None):
        df = self.df
        if bins is None:
            bins = range(0, df.read_length.max(), 100)
        mapped = df[df.reference_name != -1]
        unmapped = df[df.reference_name == -1]
        pylab.hist(mapped.read_length, bins=bins, alpha=0.5,
            label="mapped {}".format(len(mapped)), normed=True)
        pylab.hist(unmapped.read_length, bins=bins, alpha=0.5,
            label="unmapped {}".format(len(unmapped)), normed=True)
        pylab.xlabel("Isoform length")
        pylab.legend()





