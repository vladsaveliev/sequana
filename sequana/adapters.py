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
"""Utilities to manipulate adapters

Adapters removal can be performed by many different tools such as CutAdapt,
AlienTrimmer, Trimmomatic ... but they tend to use different formats from FASTA
to text files.

In this module we provide utilities to help building the input files but also
other tools to e.g. merge several list of adapters making sure there is no
duplicates.

We also store list of adapters in FASTA format, which can be read using 
the :class:`AdapterReader`::

    from sequana import sequana_data, AdapterReader
    filename = sequana_data("adapters_Nextera_PF1_220616_fwd.fa")
    ar = AdapterReader(filename)
    ar.get_adapter_by_index("N501")

"""
import pandas as pd
import pysam
from sequana.fasta import FastA
import os


def fasta_fwd_rev_to_columns(file1, file2=None, output_filename=None):
    """From two FASTA files (reverse and forward) adapters, returns 2-columns file

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
    """Convert a FASTA formatted file into adapters_ngs format

    That is a TSV file

    .. warning:: may be removed in the future.


    """
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


class Adapter():
    """Class to store one adapter

    An adapter is just a sequence from a FASTA file. It contains
    an identifier, a sequence and possibly a comment.

    .. warning:: the identifier is the FASTA identifier. It must start with a
        ">" character followed by the identifier. This is different from an
        index. 
    """
    def __init__(self, identifier, sequence, comment="nocomment"):
        self._comment = comment
        self._sequence = sequence
        self._identifier = identifier

    def _get_sequence(self):
        return self._sequence
    def _set_sequence(self, sequence):
        self._sequence = sequence
    sequence = property(_get_sequence, _set_sequence, 
        doc="R/W adapter's sequence")

    def _get_identifier(self):
        return self._identifier
    def _set_identifier(self, identifier):
        self._identifier = identifier
    identifier = property(_get_identifier, _set_identifier, 
        doc="R/W adapter's identifier")

    def _get_comment(self):
        return self._comment
    def _set_comment(self, comment):
        self._comment = comment
    comment = property(_get_comment, _set_comment, 
        doc="R/W adapter's identifier")

    def __str__(self):
        if self.comment is None:
            txt = ">%(identifier)s\n"
        else:
            txt = ">%(identifier)s\t%(comment)s\n"
        txt+= "%(sequence)s"
 
        txt = txt % {"identifier":self.identifier,
                "comment":self.comment, "sequence":self.sequence}
        return txt

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if self._comment != other.comment:
            return False
        if self._identifier != other.identifier:
            return False
        if self._sequence.upper() != other.sequence.upper():
            return False
        return True


class AdapterReader(object):
    """We use FASTA as our data structure to store adapters

    Header of the FASTA must be of the form::

        Nextera_index_N501|index_dna:N501 optional comment

    In the FASTA identifier, the part after the pipe that contains the index_dna field
    will be used to identify the index of the adapter and must be found.

    .. note:: the universal adapter does not need to have it but must be called
        *Universal_Adapter*

    .. doctest::

        >>> from sequana import sequana_data, AdapterReader
        >>> filename = sequana_data("adapters_Nextera_PF1_220616_fwd.fa")
        >>> ar = AdapterReader(filename)
        >>> candidate = ar.get_adapter_by_index("index_dna:S505")
        >>> print(candidate)
        >Nextera_index_S505|index_dna:S505
        AATGATACGGCGACCACCGAGATCTACACGTAAGGAGTCGTCGGCAGCGTC

        >>> len(ar)
        56

    """
    def __init__(self, filename):
        if isinstance(filename, str):
            # this is not large files so we load in memory all sequences/names and
            # comments once for all. This has also the adavantage that data can now
            # be changed on the fly
            fasta = FastA(filename)
            self._data = [self._to_read(this) for this in fasta]
        elif isinstance(filename, AdapterReader):
            self._data = [self._to_read(this) for this in filename._data]
        elif isinstance(filename, list):
            self._data = [self._to_read(this) for this in filename]

    def __len__(self):
        return len(self._data)

    def _get_identifiers(self):
        return [this.identifier for this in self._data]
    identifiers = property(_get_identifiers)

    def _get_seq(self):
        return [this.sequence for this in self._data]
    sequences = property(_get_seq)

    def _get_comments(self):
        return [this.comment for this in self._data]
    comments = property(_get_comments)

    def sanity_check(self):
        """Check that all identifiers are unique"""
        if len(set(self.identifiers)) != len(self.identifiers):
            import collections
            identifiers = [k for k,v in collections.Counter(self.identifiers).items() if v>1]
            msg = "Found identical identifiers in fasta sequences\n"
            msg += "Check those identifiers for duplicates %s " % identifiers
            raise ValueError(msg)

    def get_adapter_by_sequence(self, sequence):
        """Return one or several adapters that have the sequence in them

        :param str sequence: a string (ACGT letters)
        :return: Name and sequence in FASTA format that have the user *sequence*
            contained in their sequence
        """
        adapters = []
        for this in self._data:
            if sequence in this.sequence:
                this_adapter = Adapter(identifier=this.identifier, sequence=this.sequence,
                    comment=this.comment)
                adapters.append(this_adapter)
        if len(adapters) == 0:
            return None
        else:
            return adapters

    def get_adapter_by_identifier(self, text):
        """Return adapter whose identifier matches the user text

        :param index_identifier: the unique index identifier to be found. If several
            sequence do match, this is an error meaning the fasta file
            with all adapters is not correctly formatted.
        :return: the adapter that match the index_name (if any) otherwise
            returns None
        """
        adapters = []
        for this in self._data:
            if text == this.identifier or text in this.identifier.split("|"):
                this_adapter = Adapter(identifier=this.identifier, sequence=this.sequence,
                    comment=this.comment)
                adapters.append(this_adapter)

        if len(adapters) == 0:
            raise ValueError("No matching identifier found with this pattern: %s" % text)
        elif len(adapters) == 1:
            return adapters[0]
        else:
            raise ValueError("identifier %s found several times" % text)

    def get_adapter_by_index(self, index_name):
        """Return adapter corresponding to the unique index 

        :param index_name: the unique index name to be found. If several
            sequence do match, this is an error meaning the fasta file
            with all adapters is not correctly formatted.
        :return: an instance of :class:`Adapter` if index_name match an adapter, 
            returns None otherwise

        ::

            from sequana import sequana_data, AdapterReader
            filename = sequana_data("adapters_Nextera_PF1_220616_fwd.fa")
            ar = AdapterReader(filename)
            ar.get_adapter_by_identifier("index_dna:N712")

        """
        # there should be only one
        adapters = []
        for this in self._data:
            if index_name in this.identifier.split("|"):
                this_adapter = Adapter(identifier=this.identifier, sequence=this.sequence,
                    comment=this.comment)
                adapters.append(this_adapter)

        if len(adapters) == 0:
            return None
        elif len(adapters) == 1:
            return adapters[0]
        else:
            raise ValueError("index_name %s found several times" % index_name)

    def _to_read(self, this):
        from easydev import AttrDict
        d = AttrDict()
        d.sequence = this.sequence
        d.comment = this.comment
        try:
            #pysam format
            d.identifier = this.name
        except:
            # this class convention
            d.identifier = this.identifier
        return d

    def __getitem__(self, i):
        return self._data[i]

    def to_dict(self):
        d1 = [(this.identifier, [this.comment, this.sequence]) for this in self._data]
        return dict(d1)

    def __eq__(self, other):
        other = AdapterReader(other)
        d1 = self.to_dict()
        d2 = other.to_dict()
        return d1 == d2

    def reverse(self):
        """Reverse all sequences inplace

        .. doctest::

            >>> from sequana import sequana_data, AdapterReader
            >>> filename = sequana_data("adapters_Nextera_PF1_220616_fwd.fa")
            >>> filename2 = sequana_data("adapters_Nextera_PF1_220616_rev.fa")
            >>> ar = AdapterReader(filename)
            >>> ar2 = AdapterReader(filename2)
            >>> ar.reverse()
            >>> ar == ar2
            True

        """
        for this in self._data:
            this.sequence = this.sequence[::-1]

    def reverse_complement(self):
        """Reverse-complement all sequences inplace

        ::

            >>> from sequana import sequana_data, AdapterReader
            >>> filename = sequana_data("adapters_Nextera_PF1_220616_fwd.fa")
            >>> filename = sequana_data("adapters_Nextera_PF1_220616_revcomp.fa")
            >>> ar = AdapterReader(filename)
            >>> ar.reverse_complement()
            >>> ar.to_fasta()
            >>> ar == ar2

        """
        # no need for a fast implementation (less than 1000 short reads)
        from sequana.sequence import DNA
        for this in self._data:
            this.sequence = DNA(this.sequence[:]).get_reverse_complement()

    def to_fasta(self, filename):
        """Save sequences into fasta file

        .. todo:: if no comment, currently None is written. Would be better to
            simply skip it
        """
        with open(filename, "w") as fout:
            for i, this in enumerate(self._data):
                if i>0:
                    fout.write("\n")
                fout.write(">%s\t%s\n%s" % (this.identifier, this.comment, this.sequence))


class FindAdaptersFromIndex(object):
    """Used to identify one or two indices from an IndexMapper structure

    An IndexMapper structure is just a CSV file that stores
    the mapping between a sample name and one or two indices.

    The columns are project, index1, index2 and sample_name. 

    Used by sequana main script top build the adapter files for multi-sample
    projects. 

    """
    def __init__(self, index_mapper, adapters, mode="revcomp"):
        """.. rubric:: Constructor

        :param str index_mapper: filename of a CSV file that has the following
            header projet,index1,index2,sample_name (index2 is optional)

        sample name is the string that starts the filename. From
        the sample name, get the indices 1 and 2 and then the corresponding
        adapter.

        """
        if mode != "revcomp":
            raise ValueError("only mode revcomp implemented. ")
        #self.index_mapper = pd.read_csv(index_mapper, delim_whitespace=True)[columns_in]
        try:
            columns_in = ['sample_name', 'index1', 'index2']
            self.index_mapper = pd.read_csv(index_mapper, sep=",", dtype=str)[columns_in]
            self.mode = 'multi_index'
        except:
            columns_in = ['sample_name', 'index1']
            self.index_mapper = pd.read_csv(index_mapper, sep=",", dtype=str)[columns_in]
            self.mode = 'single_index'

        self.index_mapper.columns = columns_in
        self.index_mapper.set_index('sample_name', inplace=True)

        self.adapters = adapters
        if adapters == "Nextera":
            from sequana import sequana_data
            file1 = sequana_data("adapters_Nextera_PF1_220616_fwd.fa", 
                "data/adapters")
            file2 = sequana_data("adapters_Nextera_PF1_220616_revcomp.fa",
                "data/adapters")
        elif adapters == "PCRFree":
            from sequana import sequana_data
            file1 = sequana_data("adapters_PCR-free_PF1_220616_fwd.fa", 
                "data/adapters")
            file2 = sequana_data("adapters_PCR-free_PF1_220616_revcomp.fa",
                "data/adapters")
        else:
            raise ValueError("Unknown set of adapters. Use Nextera or PCRFree (note the small/big caps)")

        self._adapters_fwd = AdapterReader(file1)
        self._adapters_rev = AdapterReader(file2)  # !!! revcomp

    def _get_samples(self):
        return list(self.index_mapper.index)
    sample_names = property(_get_samples)

    def get_indices(self, sample_name):
        if sample_name not in self.index_mapper.index:
            raise ValueError("%s not valid. Use one of %s" % (sample_name,
                                                              self.sample_names))
        return self.index_mapper.ix[sample_name]

    def get_adapters(self, sample_name, include_universal=True,
            include_transposase=True):

        indices = self.get_indices(sample_name)

        if self.mode == "single_index" or self.mode == "multi_index":
            res = {'index1': {}}
            index1 = "index_dna:" + indices.ix['index1']
            #index1 = indices.ix['index1']
            res['index1']['fwd'] = self._adapters_fwd.get_adapter_by_index(index1)
            res['index1']['rev'] = self._adapters_rev.get_adapter_by_index(index1)

        if self.mode == "multi_index":
            res['index2'] = {}
            index2 = "index_dna:" + str(indices.ix['index2'])
            #index2 =  indices.ix['index2']
            res['index2']['fwd'] = self._adapters_fwd.get_adapter_by_index(index2)
            res['index2']['rev'] = self._adapters_rev.get_adapter_by_index(index2)

        if include_universal:
            res['universal'] = {}
            res['universal']['fwd'] = self._adapters_fwd.get_adapter_by_identifier(
                'Universal_Adapter')
            res['universal']['rev'] = self._adapters_rev.get_adapter_by_identifier(
                'Universal_Adapter')

        if include_transposase and self.adapters == "Nextera":
            res['transposase'] = {}
            res['transposase']['fwd'] = str(self._adapters_fwd.get_adapter_by_identifier(
                'Nextera_transposase_seq_1'))
            res['transposase']['fwd'] += "\n" + str(self._adapters_fwd.get_adapter_by_identifier(
                'Nextera_transposase_seq_2'))
            res['transposase']['rev'] = str(self._adapters_rev.get_adapter_by_identifier(
                'Nextera_transposase_seq_1'))
            res['transposase']['rev'] += "\n"+str(self._adapters_rev.get_adapter_by_identifier(
                'Nextera_transposase_seq_1'))


        return res

    def save_adapters_to_fasta(self, sample_name, include_universal=True,
            include_transposase=True, output_dir='.'):
        """Get index1, index2 and uiversal adapter"""
        adapters = self.get_adapters(sample_name, include_universal=include_universal)

        file_fwd = output_dir + os.sep + "%s_adapters_fwd.fa"% sample_name
        with open(file_fwd, "w") as fout:
            if include_transposase and self.adapters == "Nextera":
                fout.write(str(adapters['transposase']['fwd'])+"\n")

            if include_universal:
                fout.write(str(adapters['universal']['fwd'])+"\n")

            fout.write(str(adapters['index1']['fwd'])+"\n")
            if self.mode == "multi_index":
                fout.write(str(adapters['index2']['fwd'])+"\n")

        file_rev = output_dir + os.sep + "%s_adapters_revcomp.fa" % sample_name
        with open(file_rev, "w") as fout:
            if include_transposase and self.adapters == "Nextera":
                fout.write(str(adapters['transposase']['rev'])+"\n")
            if include_universal:
                fout.write(str(adapters['universal']['rev'])+"\n")
            fout.write(str(adapters['index1']['rev'])+"\n")
            if self.mode == "multi_index":
                fout.write(str(adapters['index2']['rev'])+"\n")

        return file_fwd, file_rev


class AdapterRemoval(object):
    """Possible data structure for future usage. Not used yet"""
    def __init__(self, file1, file2=None, adapter1=None, adapter2=None,
        quality_cutoff1=30, quality_cutoff2=None):

        self.file1 = file1
        self.file2 = file2
        self.qual1 = quality_cutoff1
        self.qual2 = quality_cutoff2
        self.adapter1 = adapter1
        self.adapter2 = adapter2

    def save_adapters_to_fasta(self, sample_name, include_universal=True, output_dir='.'):
        """Get index1, index2 and uiversal adapter"""
        adapters = self.get_adapters(sample_name, include_universal=include_universal)

        file_fwd = output_dir + os.sep + "%s_adapters_fwd.fa"% sample_name
        with open(file_fwd, "w") as fout:
            if include_universal:
                fout.write(str(adapters['universal']['fwd'])+"\n")
            fout.write(str(adapters['index1']['fwd'])+"\n")
            if self.mode == "multi_index":
                fout.write(str(adapters['index2']['fwd'])+"\n")

        file_rev = output_dir + os.sep + "%s_adapters_rev.fa" % sample_name
        with open(file_rev, "w") as fout:
            if include_universal:
                fout.write(str(adapters['universal']['rev'])+"\n")
            fout.write(str(adapters['index1']['rev'])+"\n")
            if self.mode == "multi_index":
                fout.write(str(adapters['index2']['rev'])+"\n")

        return file_fwd, file_rev


class AdapterRemoval(object):
    """Possible data structure for future usage. Not used yet"""
    def __init__(self, file1, file2=None, adapter1=None, adapter2=None,
        quality_cutoff1=30, quality_cutoff2=None):

        self.file1 = file1
        self.file2 = file2
        self.qual1 = quality_cutoff1
        self.qual2 = quality_cutoff2
        self.adapter1 = adapter1
        self.adapter2 = adapter2
