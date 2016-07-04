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
from sequana.fasta import FastA
import os

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
            self.load_fasta(sequana_data("%s" % v, "data"))

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

class AdapterReader(FastA):
    """We use FastA as our data structure to store adapters

    Header of the FASTA must be of the form::

        Nextera_index_N501|index_dna:N501

    with a *|index_dna:* string followed by the index tag

    .. note:: the universal adapter does not need to have it but is called
        Universal_Adapter


    .. doctest::

        >>> from sequana import sequana_data, AdapterReader
        >>> ar = AdapterReader(sequana_data("adapters_Nextera_PF1_220616_fwd.fa"))
        >>> candidate = ar.get_get_adapter_by_index("S505")
        >>> print(candidate)
        '>Nextera_index_S505|index_dna:S505
        AATGATACGGCGACCACCGAGATCTACACGTAAGGAGTCGTCGGCAGCGTC'


    """
    def __init__(self, filename):
        super(AdapterReader, self).__init__(filename)
        # TODO check uniqueness of the indices

    def get_adapter_by_sequence(self, sequence):
        """Return one or several adapters that have the sequence in them

        :param str sequence: a string (ACGT letters)
        :return: Name and sequence in FASTA format that have the user *sequence*
            contained in their sequence
        """
        adapters = [str(this) for this in self if sequence in this.sequence]
        if len(adapters) == 0:
            return None
        else:
            return "\n".join(adapters)

    def get_adapter_by_name(self, text):
        """Return adapter whose name matches the user text

        :param index_name: the unique index name to be found. If several
            sequence do match, this is an error meaning the fasta file
            with all adapters is not correctly formatted.
        :return: the adapter that match the index_name (if any) otherwise
            returns None
        """
        adapters = [str(this) for this in self if text in this.name]
        if len(adapters) == 0:
            return None
        elif len(adapters) == 1:
            return adapters[0]
        else:
            raise ValueError("name %s found several times" % text)

    def get_adapter_by_index(self, index_name):
        """

        :param index_name: the unique index name to be found. If several
            sequence do match, this is an error meaning the fasta file
            with all adapters is not correctly formatted.
        :return: the adapter that match the index_name (if any) otherwise
            returns None
        """
        """Return FASTA corresponding to the index"""
        # there should be only one
        adapters = [str(this) for this in self if index_name in this.name]
        if len(adapters) == 0:
            return None
        elif len(adapters) == 1:
            return adapters[0]
        else:
            raise ValueError("index_name %s found several times" % index_name)


class FindAdaptersFromIndex(object):
    """

    """
    def __init__(self, index_mapper, adapters='Nextera'):
        """.. rubric:: Constructor

        :param str index_mapper: filename of a CSV file that has the following
            header Projet,Index1,Index2,sample name,Index1,Index2

        sample name is the string that starts the filename. From
        the sample name, get the indices 1 and 2 ann then the corresponding
        adapters.

        """
        self.index_mapper = pd.read_csv(index_mapper)[["sample name", 'Index1.1', 'Index2.1']]
        self.index_mapper.columns = ["sample", "index1", "index2"]
        self.index_mapper.set_index('sample', inplace=True)

        if adapters == "Nextera":
            from sequana import sequana_data
            file1 = sequana_data("adapters_Nextera_PF1_220616_fwd.fa", "data")
            file2 = sequana_data("adapters_Nextera_PF1_220616_rev.fa", "data")

            self._adapters_fwd = AdapterReader(file1)
            self._adapters_rev = AdapterReader(file2)
        else:
            raise NotImplementedError

    def _get_samples(self):
        return list(self.index_mapper.index)
    sample_names = property(_get_samples)

    def get_indices(self, sample_name):
        if sample_name not in self.index_mapper.index:
            raise ValueError("%s not valid. Use one of %s" % (sample_name,
                                                              self.sample_names))
        return self.index_mapper.ix[sample_name]

    def get_adapters(self, sample_name, include_universal=True):
        indices = self.get_indices(sample_name)

        res = {'index1': {}, 'index2': {}}

        index1 = indices.ix['index1']
        res['index1']['fwd'] = self._adapters_fwd.get_adapter_by_index(index1)
        res['index1']['rev'] = self._adapters_rev.get_adapter_by_index(index1)


        index2 = indices.ix['index2']
        res['index2']['fwd'] = self._adapters_fwd.get_adapter_by_index(index2)
        res['index2']['rev'] = self._adapters_rev.get_adapter_by_index(index2)

        if include_universal:
            res['universal'] = {}
            res['universal']['fwd'] = self._adapters_fwd.get_adapter_by_name(
                'Universal_Adapter')
            res['universal']['rev'] = self._adapters_rev.get_adapter_by_name(
                'Universal_Adapter')
        return res

    def save_adapters_to_csv(self, sample_name, include_universal=True, output_dir='.'):
        """Get index1, index2 and uiversal adapter"""
        adapters = self.get_adapters(sample_name, include_universal=include_universal)

        file_fwd = output_dir + os.sep + "%s_adapters_fwd.fa"% sample_name
        with open(file_fwd, "w") as fout:
            if include_universal:
                fout.write(adapters['universal']['fwd']+"\n")
            fout.write(adapters['index1']['fwd']+"\n")
            fout.write(adapters['index2']['fwd']+"\n")

        file_rev = output_dir + os.sep + "%s_adapters_rev.fa" % sample_name
        with open(file_rev, "w") as fout:
            if include_universal:
                fout.write(adapters['universal']['rev']+"\n")
            fout.write(adapters['index1']['rev']+"\n")
            fout.write(adapters['index2']['rev']+"\n")

        return file_fwd, file_rev


class AdapterRemoval(object):
    def __init__(self, file1, file2=None, adapter1=None, adapter2=None,
        quality_cutoff1=30, quality_cutoff2=None):

        self.file1 = file1
        self.file2 = file2
        self.qual1 = quality_cutoff1
        self.qual2 = quality_cutoff2
        self.adapter1 = adapter1
        self.adapter2 = adapter2
