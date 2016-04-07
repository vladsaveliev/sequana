"""




"""
import os
import pandas as pd
import pylab
import pysam
from reports import HTMLTable

"""
#http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42#
#http://biofinysics.blogspot.fr/2014/05/how-does-bowtie2-assign-mapq-scores.html
#https://gitlab.univ-nantes.fr/a-slide/ContaVect/blob/9a411abfa720064c205c5f6c811afdfea206ed12/pyDNA/pySamTools/Bam.py
# pysamtools
# pysam uses htlib behing the scene and is very fast
# pysam works great for BAM file but with SAM, it needs to read the file after
# each compete iteration, which is not very useful

Interesting commands::

    samtools flagstat contaminant.bam
    samtools stats contaminant.bam
"""

class BAM(pysam.AlignmentFile):
    """

    mode rb for bam files

    .. doctest::

        >>> from sequana import BAM, sequana_data
        >>> b = BAM(sequana_data("test.bam"))
        >>> len(b)
        1000

    """
    def __init__(self, filename, mode="rb", *args):
        # works for py27 but not py35 probably a missing __init__  or __new__ in
        # AlignmentFile class. See e.g., stackoverflow/questions/
        # 26653401/typeerror-object-takes-no-parameters-but-only-in-python-3
        #super(BAM, self).__init__(filename, mode, *args)
        pysam.AlignmentFile.__init__(filename, mode, *args)

    def get_read_names(self):
        """Return the reads' names"""
        self.reset()
        names = [this.qname for this in self]
        self.reset()
        return names

    def iter_unmapped_reads(self):
        self.reset()
        unmapped = (this.qname for this in self if this.is_unmapped)
        self.reset()
        return unmapped

    def iter_mapped_reads(self):
        self.reset()
        mapped = (this.qname for this in self if this.is_unmapped is False)
        self.reset()
        return mapped

    def __len__(self):
        self.reset()
        N = len([x for x in self])
        self.reset()
        return N

    def get_stats(self):
        """Return basic stats about the reads



        """

        """#See samtools stats

        # 1526795 + 0 in total (QC-passed reads + QC-failed reads)
        13 + 0 secondary
        0 + 0 supplementary
        0 + 0 duplicates
        3010 + 0 mapped (0.20% : N/A)
        1526782 + 0 paired in sequencing
        763391 + 0 read1
        763391 + 0 read2
        2700 + 0 properly paired (0.18% : N/A)
        2976 + 0 with itself and mate mapped
        21 + 0 singletons (0.00% : N/A)
        0 + 0 with mate mapped to a different chr
        0 + 0 with mate mapped to a different chr (mapQ>=5)
        """

        d = {}
        d['total_reads'] = len(list(self.iter_unmapped_reads()))
        d['mapped_reads'] = len(list(self.iter_mapped_reads()))
        d['unmapped_reads'] = len(list(self.iter_unmapped_reads()))
        d['contamination [%]'] = float(d['mapped_reads']) /float(d['unmapped_reads']) 
        d['contamination [%]'] *= 100
        return d

    def get_flags(self):
        self.reset()
        flags = [s.flag for s in self]
        self.reset()
        return flags

    def get_flags_as_df(self):
        """

        .. doctest::

            >>> from sequana import BAM, sequana_data
            >>> b = BAM(sequana_data('test.bam'))
            >>> df = b.get_flags_as_df()
            >>> df.sum()
            1       1000
            2        484
            4          2
            8          2
            16       499
            32       500
            64       477
            128      523
            256       64
            512        0
            1024       0
            2048       0
            dtype: int64

        """
        flags = self.get_flags()
        data = [(this, [flag&this for flag in flags]) 
            for this in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]]
        import pandas as pd
        df = pd.DataFrame(dict(data))
        df = df > 0
        return df

    def plot_bar_flags(self, logy=True, fontsize=16, filename=None):
        """

        .. plot::
            :include-source:

            from sequana import BAM, sequana_data
            b = BAM(sequana_data('test.bam'))
            b.plot_bar_flags()


        """
        df = self.get_flags_as_df()
        df = df.sum()
        pylab.clf()
        if logy is True:
            df.plot(kind='bar', logy=logy, grid=True)
        else:
            df.plot(kind='bar', grid=True)
        pylab.xlabel("flags", fontsize=fontsize)
        pylab.ylabel("count", fontsize=fontsize)
        try:
            pylab.tight_layout()
        except:
            pass
        if filename:
            pylab.savefig(filename)

    def to_fastq(self):
        """

        .. todo::

        """
        raise NotImplementedError

    def to_sam(self):
        """

               .. todo::

        """
        raise NotImplementedError

    def head(self, N=100):
        """

               .. todo::

        """
        # Export the first top alignment into a new BAM file ?
        raise NotImplementedError

    def get_mapq_as_df(self):
        """Return dataframe with mapq for each read


        """
        self.reset()
        df = pd.DataFrame({'mapq': [this.mapq for this in self]})
        self.reset()
        return df

    def plot_bar_mapq(self, fontsize=16, filename=None):
        """

            .. plot::
                :include-source:

                from sequana import BAM, sequana_data
                b = BAM(sequana_data('test.bam'))
                b.plot_bar_mapq()

            """

        df = self.get_mapq_as_df()
        df.plot(kind='hist', bins=60, legend=False, grid=True, logy=True)
        pylab.xlabel("MAPQ", fontsize=fontsize)
        try:
            pylab.tight_layout()
        except:
            pass
        if filename:
            pylab.savefig(filename)



class Alignment(object):
    """Helper class to retrieve info about Alignment

    Takes an alignment as read by :class:`BAM` and provide simplified version
    of pysam.Alignmennt class.

    ::

        >>> from sequana.bamtools import Alignment
        >>> from sequana import BAM, sequana_data
        >>> b = BAM(sequana_data("test.bam"))
        >>> segment = next(b)
        >>> align = Alignment(segment)
        >>> align.as_dict()
        >>> align.FLAG
        353



    The original data is stored in hidden attribute :attr:`_data` and the
    following values are available as attributes or dictionary:


    * QNAME: a query template name. Reads/segment having same QNAME come from the
      same template. A QNAME set to `*` indicates the information is unavailable.
      In a sam file, a read may occupy multiple alignment 
    * FLAG: combination of bitwise flags. See :class:`SAMFlags`
    * RNAME: reference sequence
    * POS
    * MAPQ: mapping quality if segment is mapped. equals -10 log10 Pr
    * CIGAR: 
    * RNEXT: reference sequence name of the primary alignment of the NEXT read
      in the template
    * PNEXT: position of primary alignment
    * TLEN: signed observed template length
    * SEQ: segment sequence
    * QUAL: ascii of base quality

    .. note:: A segment is a contiguous sequence. A read is a sequence that
        may consist of multiple segments.


    """
    def __init__(self, alignment):
        """.. rubric:: constructor

        :param alignment: alignment instance from :class:`BAM`


        """
        self._data = alignment
        d = self.as_dict()
        for key in d.keys():
            setattr(self, key, d[key])

    def as_dict(self):
        d = {}
        s = self._data
        d['QNAME'] = s.qname
        d['FLAG'] = s.flag
        d['RNAME'] = s.rname
        d['POS'] = s.pos
        d['MAPQ'] = s.mapq
        d['CIGAR'] = s.cigar
        d['PNEXT'] = s.pnext
        d['RNEXT'] = s.rnext
        d['TLEN'] = s.tlen
        d['SEQ'] = s.seq
        d['QUAL'] = s.qual
        return d




class SAMFlags(object):
    """

    .. doctest::

        >>> from sequana import SAMFlags
        >>> sf = SAMFlags(257)
        >>> sf.get_flags()
        [1, 256]


    You can also print the bits and their description::

        print(sf)

    ======= ====================================================================
    bit     Meaning/description
    ======= ====================================================================
    1       template having multiple segments in sequencing
    2       each segment properly aligned according to the aligner
    4       segment unmapped
    8       next segment in the template unmapped
    16      SEQ being reverse complemented
    32      SEQ of the next segment in the template being reverse complemented
    64      the first segment in the template
    128     the last segment in the template
    256     secondary alignment
    512     not passing filters, such as platform/vendor quality controls
    1024    PCR or optical duplicate
    2048    supplementary alignme
    ======= ====================================================================

    """
    def __init__(self, value=4095):
        self.value = value
        self._flags = {
            1: "template having multiple segments in sequencing",
            2: "each segment properly aligned according to the aligner",
            4: "segment unmapped",
            8: "next segment in the template unmapped",
            16: "SEQ being reverse complemented",
            32: "SEQ of the next segment in the template being reverse complemented", 
            64: "the first segment in the template",
            128: "the last segment in the template",
            256: "secondary alignment",
            512: "not passing filters, such as platform/vendor quality controls",
            1024: "PCR or optical duplicate",
            2048: "supplementary alignme"}

    def get_meaning(self):
        """

        Return the meaning/description (sorted by bit i.e. 1,2,4,8...2048)
        """
        return [self._flags[k] for k in sorted(self._flags.keys())]

    def get_flags(self):
        """Return the individual bits included in the flag"""
        flags = []
        for this in sorted(self._flags.keys()):
            if self.value & this:
                flags.append(this)
        return flags

    def __str__(self):
        txt = ""
        for this in sorted(self._flags.keys()):
            if self.value & this:
                txt += "%s: %s\n" % (this, self._flags[this])
        return txt


from .report_main import BaseReport
class BAMReport(BaseReport):
    """

    ::

        from sequana import BAM, sequana_data, BAMReport
        b = BAM(sequana_data("test.bam"))

        r = BAMReport()
        r.set_data(b)
        r.create_report()

        # report/bam.html is now available

    """
    def __init__(self, jinja_template="bam", output_filename="bam.html",
                 directory="report", **kargs):
        super(BAMReport, self).__init__(jinja_template, output_filename,
            directory, **kargs)

        self.jinja['title'] = "Bam Report"

    def set_data(self, data):
        self.bam = data

    def parse(self):
        self.jinja['alignment_count'] = len(self.bam)

        # first, we store the flags
        df = self.bam.get_flags_as_df().sum()
        df = df.to_frame()
        df.columns = ['counter']
        sf = SAMFlags()
        df['meaning'] = sf.get_meaning()
        df = df[['meaning', 'counter']]
        html = HTMLTable(df).to_html(index=True)
        self.jinja['flags_table'] = html

        # create the bar plot with flags
        self.bam.plot_bar_flags(logy=True, filename=self.directory + os.sep +
                                                    "bar_flags_logy.png")
        self.bam.plot_bar_flags(logy=False, filename=self.directory + os.sep +
                                                     "bar_flags.png")

        self.bam.plot_bar_mapq(filename=self.directory + os.sep + "bar_mapq.png")


