import re
import string
from collections import Counter

from sequana.fasta import FastA


class Sequence(object):
    """Abstract base classe for other specialised sequences such as DNA.


    Sequenced is the base class for other classes such as :class:`DNA` and
    :class:`RNA`.

    ::

        from sequana import Sequence
        s = Sequence("ACGT")
        s.stats()
        s.get_complement()

    .. note:: You may use a Fasta file as input (see constructor)


    """
    def __init__(self, sequence, complement_in=b"ACGT", complement_out=b"TGCA",
                 letters="ACGT"):
        """.. rubric:: Constructor

        A sequence is just a string stored in the :attr:`sequence` attribute. It
        has properties related to the type of alphabet authorised.

        :param str sequence: May be a string of a Fasta File, in which case only
            the first sequence is used.
        :param complement_in:
        :param complement_out:
        :param letters: authorise letters. Used in :meth:`check` only.

        .. todo:: use counter only once as a property

        """
        if sequence.endswith(".fa") or sequence.endswith(".fasta"):
            fasta = FastA(sequence)
            sequence = fasta.next().sequence
        else: # assume correct string sequence
            pass

        self._data = sequence
        try:
            self._translate = string.maketrans(complement_in, complement_out)
        except:
            self._translate = bytes.maketrans(complement_in, complement_out)
        self._letters = letters

    def _get_sequence(self):
        return self._data
    sequence = property(_get_sequence)

    def get_complement(self):
        """Return complement """
        return self._data.translate(self._translate)

    def get_reverse_complement(self):
        """Return reverse complement """
        return self.get_complement()[::-1]

    def get_reverse(self):
        """Return reverse sequence"""
        return self._data[::-1]

    def complement(self):
        """Alias to :meth:`get_complement`"""
        self._data = self.get_complement()

    def reverse(self):
        """Alias to :meth:`get_reverse`"""
        self._data = self.get_reverse()

    def reverse_complement(self):
        """Alias to get_reverse_complement"""
        self._data = self.get_reverse_complement()

    def check(self):
        """Check that all letters are valid"""
        counter = Counter(self._data).keys()
        for key in counter:
            if key not in self._letters:
                raise ValueError("Found unexpected letter in the sequence (%s)" % key)

    def __len__(self):
        return len(self._data)

    def gc_content(self):
        """Return mean GC content"""
        c = Counter(self._data)
        ratio = (c['G'] + c['C']) / len(self.sequence)
        return ratio

    def stats(self):
        """Return basic stats about the number of letters"""
        from collections import Counter
        return Counter(self.sequence)

    def get_occurences(self, pattern, overlap=False):
        """Return position of the input pattern in the sequence

        ::

            >>> from sequana import Sequence
            >>> s = Sequence('ACGTTTTACGT')
            >>> s.get_occurences("ACGT")
            [0, 7]

        """
        if overlap is False:
            res = [m.start() for m in re.finditer(pattern, self.sequence)]
        elif overlap is True:
            res = [m.start() for m in re.finditer('(?=%s)'%pattern, self.sequence)]
        return res

        # reverse find-all without overlaps, you can combine positive and
        # negative lookahead into an expression like this:
        #res = [m.start() for m in re.finditer('(?=%s)(?!.{1,%d}%s)' % (search, 
        #    len(pattern)-1, pattern), 'ttt')]


class DNA(Sequence):
    """Simple DNA class


        >>> d = DNA("ACGTTTT")
        >>> d.complement
        >>> d.reverse_complement

    """
    def __init__(self, sequence):
        super(DNA, self).__init__(sequence, complement_in=b"ACGTacgt",
            complement_out=b"TGCAtgca", letters="ACGTacgtNn")


class RNA(Sequence):
    """Simple RNA class


        >>> d = RNA("ACGUUUU")
        >>> d.complement
        >>> d.reverse_complement

    """
    def __init__(self, sequence):
        super(RNA, self).__init__(sequence, complement_in=b"ACGUacgu",
            complement_out=b"UGCAugca", letters="ACGUacguNn")
