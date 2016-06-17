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
"""Adapters

Adapters removal can be performed by many different tools such as CutAdapt,
AlienTrimmer, Trimmomatic ... but they tend to use different formats from FASTA
to text files. Besides, list of adapters are usually in FASTA formats. 

In this module we provide utilities to help building the input files but also
other tools to e.g. merge several list of adapters making sure there is no
duplicates.

This module may also be used (WIP) by the kraken module to provide a list
of adapters in FASTA format but with specific annotation. Consider for instance
this adapter in FASTA format::

    >NextFlex_PCR_Free_adapter1   NextFlex_PCR_Free_adapter1
    GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG

The name is not standard here. If you want use this FASTA in
other tools, it may fail so we can added the missing bits that is::

    >NextFlex_PCR_Free_adapter1|kraken:taxid|10000001   NextFlex_PCR_Free_adapter1
    GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG

"""
import pandas as pd
import pysam


def fasta_fwd_rev_to_columns(file1, file2=None, output_filename=None):
    """Convert reverse and forward adapters (FASTA format) into 2-columns file

    This is useful for some tools related to adapter removal that takes as input
    this kind of fomat

    :param str filename1: FASTA format
    :param str filename2: FASTA format (optional)

    The files must have a one-to-one mapping
    """
    f1 = pysam.FastxFile(file1)
    if output_filename is not None:
        fout = open(output_filename, "w")
    if file2:
        f2 = pysam.FastxFile(file2)
        for read1, read2 in zip(f1, f2):
            txt = "%s %s" % (read1.sequence, read2.sequence)
            if output_filename is None:
                print(txt)
            else:
                fout.write(txt+"\n")
    else:
        for read1 in f1:
            txt = "%s" % read1.sequence
            if output_filename is None:
                print(read1.sequence)
            else:
                fout.write(txt+"\n")
    if output_filename is not None:
        fout.close()


def adapters_to_clean_ngs(input_filename, output_filename="adapters_ngs.txt"):
    with open(input_filename, 'r') as fh1:
        data1 = fh1.readlines()

    count = 0
    with open(output_filename, "w") as fout:
        for line in data1:
            line = line.strip().strip("\n")
            if line.startswith('>'):
                pass
            else:
                data = "adapter_%s\t%s\t0.5\t31\t10\t0\t0\n"% (count+1, line)
                fout.write(data)
                count+=1


def adapter_removal_parser(filename):
    """Parses output of AdapterRemoval


    .. doctest::

        >>> from sequana import adapters, sequana_data
        >>> data = sequana_data("test_adapter_removal_output.txt", "testing")
        >>> results = adapters.adapter_removal_parser(data)
        >>> results["adapter1"]
        'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'


    """
    results = {}

    with open(filename, "r") as fin:
        lines = fin.readlines()
        for line in lines:
            if line.startswith("  --adapter"):
                lhs, rhs = line.split(":")
                name = lhs.strip().replace("-", "")
                sequence = rhs.strip()
                results[name] = sequence
    return results


class AdapterDB(object):
    """Utility used in Kraken pipeline

    The name of the Fasta should be formatted as::

        >Name|kraken:taxid|id

    where id is a number starting with 1000000

    .. warning:: this convention is adopted for now but may change.

    """
    def __init__(self, filename=None):

        self.df = pd.DataFrame(columns=["name", "sequence",
            "comment", "identifier", "filename"])

        if filename:
            self.load_fasta(filename)

    def load_all(self):
        from sequana.resources.data import adapters as dict_adapters
        from sequana import sequana_data
        for k,v in dict_adapters.items():
            self.load_fasta(sequana_data("data/%s" % v))

    def load_fasta(self, filename):
        from sequana.fasta import FastA
        adapters = FastA(filename)
        self.records = []
        for adapter in adapters:
            identifier = adapter.name.split("|")[2]
            record = {
                'name': adapter.name,
                "sequence": adapter.sequence,
                "comment":adapter.comment,
                "identifier":identifier,
                "filename":filename}
            self.records.append(record)

        self.df = self.df.append(self.records)
        self.df.reset_index(drop=True, inplace=True)

        # check that identifiers should be unique
        if len(self.df) > len(self.df.identifier.unique()):
            print("Warn: there are duplicated identifiers in the adapters")

    def get_name(self, identifier):
        name =  self.df[self.df.identifier == str(identifier)].comment
        if len(name) == 1:
            name = list(name)[0]
        return name






