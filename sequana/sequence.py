import string


class Sequence(object):
    def __init__(self, sequence, complement_in="ACGT", complement_out="TGCA", letters=""):
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
        return self._data.translate(self._translate)

    def get_reverse_complement(self):
        return self.get_complement()[::-1]

    def get_reverse(self):
        return self._data[::-1]

    def complement(self):
        self._data = self.get_complement()

    def reverse(self):
        self._data = self.get_reverse()

    def reverse_complement(self):
        self._data = self.get_reverse_complement()

    def check(self):
        import collections
        counter = collections.Counter(self._data).keys()
        for key in counter:
            if key not in self._letters:
                raise ValueError("Found unexpected letter in the sequence (%s)" % key)

class DNA(Sequence):
    """Simple DNA class


        >>> d = DNA("ACGTTTT")
        >>> d.complement
        >>> d.reverse_complement

    """
    def __init__(self, sequence):
        super(DNA, self).__init__(sequence, complement_in=b"ACGTacgt",
            complement_out=b"TGCAtgca", letters="ACGTacgtNn")

