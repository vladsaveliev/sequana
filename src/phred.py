"""

phred scales indicates the range of of characters.
characters goes from ! to ~ that is tfrom 33 to 126in ascii tqble. Characters before ! caused trouble (e.g. white spaces). This scale is the Sanger scale. 2 other scales could be used ranging from 59 to 126 (illumina 1) and from 64 to 126 (illumina 1.3+).


http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/

Here qre the offset to use::

Sanger : 33       range 0 to 93
Solexa: 64        range -5 to 62
illumina1.3+: 64  range 0 to 62


::

    p = linspace(0,1,1000)
    plot(p, [phred.phred_quality_score_sanger(this)
        for this in p], color='red' )
    plot(p, [phred.phred_score_solexa(this)
        for this in p], color='black' )



"""
import numpy as np
import pylab

quality = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP"""
quality += """QRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""


from math import log10


__all__ = ['Quality']


def proba_to_quality_sanger(pe):
    """A value beitween 0 and 93

pe is the probability of error.
Q is the quality score.

a high probability of error (0.99) gives Q=0
q low proba of errors (0.05) gives Q = 13
q low proba of errors (0.01) gives Q = 20

    """
    if pe > 1:
        pe = 1
    if pe < 1e-90:
        pe = 1e-90
    Qs = -10 * log10(pe)
    if Qs > 93:
        Qs = 93
    return Qs

def quality_to_proba_sanger(quality):
    return 10**(quality/-10.)


def proba_to_quality_solexa(pe):
    """prior v1.3 (ref: wikipedia
    https://en.wikipedia.org/wiki/FASTQ_format
    """
    if pe > 1:
        pe = 1
    if pe <1e-90:
        pe = 1e-90
    Qs = -10 * log10(pe/(1-pe))
    if Qs > 62:
        Qs = 62
    if Q < -5:
        Q = -5
    return Qs


def quality_to_proba_sanger(quality):
    return 10**(quality/-10.)


def quality_solexa_to_quality_sanger(qual):
    return 10 * log10(10**(qual/10.) + 1 )


def quality_sanger_to_quality_solexa(qual):
    return 10 * log10(10**(qual/10.) - 1 )


def ascii_to_quality(character, phred=33):
    """

    ::

        >>> ascii_to_quality("!")
        0
        >>> ascii_to_quality("~")
        93
    """
    return ord(character) - phred


def quality_to_ascii(quality, phred=33):
    """
    ::

         >>> quality_to_ascii(65)
         b

    """
    return chr(quality+33)


class Quality(object):
    """


    .. plot::

        >>> from sequana.phread import Quality
        >>> q = Quality('BCCFFFFFHHHHHIIJJJJJJIIJJJJJJJJFH')
        >>> q.plot()
        >>> q.mean_quality()
        35
        


    """
    def __init__(self, seq, offset=33):
        self.seq = seq
        self.offset = offset

    def _get_quality(self):
        # bytearray conversion is required since ord() function
        # does not handle bytes from py3
        return [x - self.offset for x in bytearray(self.seq)]
    quality = property(_get_quality)

    def _get_mean_quality(self):
        return np.mean(self.quality)
    mean_quality = property(_get_mean_quality)

    def plot(self):
        pylab.plot(self.quality)
        pylab.xlabel('base position')
        pylab.ylabel('Quality per base')
        pylab.grid()



# this should qlso be correct for Illumina 1.8+
class QualitySanger(Quality):
    def __init__(self, seq, offset=33):
        super(QualitySanger, self).__init__(seq, offset)


class QualitySolexa(Quality):
    def __init__(self, seq, offset=64):
        super(QualitySanger, self).__init__(seq, offset)


















