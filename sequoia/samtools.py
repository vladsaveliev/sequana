"""



http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42

http://biofinysics.blogspot.fr/2014/05/how-does-bowtie2-assign-mapq-scores.html

"""


import pandas as pd
import pylab


class SAMTools(object):
    """

    Header of a SAM file has N lines starting with '@' character.
    Using comment='@' and header None does not make the job.Using skiprows=N
    works but requires to know the number of lines in the header.


    .. todo:: read by chunk size for large files ?
    .. todo:: read two files for comparison ?
    """
    def __init__(self, filename):
        self.filename = filename
        skiprows = self._guess_header_length()
        print('Found a header of %s lines. skipped' % skiprows)
        self.df = pd.read_csv(filename, sep="\t", header=None, 
                skiprows=skiprows)

    def _guess_header_length(self):
        fh = open(self.filename, 'r')
        N = 0
        skiprows = 0
        while N < 1000:
            line = fh.readline()
            if line.startswith('@'):
                skiprows += 1
            else:
                N = 1000
            N += 1
        return skiprows

    def get_mapq(self):
        data = self.df[4]
        return data

    def plot_mapq_distribution(self):
        """Plot distribution of MAPQ scores (fifth column of SAM)
        
        
        The maximum MAPQ value that Bowtie 2 generates is 42 (though it doesn't
        say this anywhere in the documentation). In contrast, the maximum MAPQ
        value that BWA will generate is 37 (though once again, you -
        frustratingly - won't find this information in the manual).

        :reference: http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42
        """
        data = self.get_mapq()
        Nmax = max(data)
        pylab.hist(data, bins=Nmax, normed=True)
        pylab.grid(True)
        pylab.xlim([0, Nmax])
        pylab.xlabel('MAPQ score')
        pylab.ylabel('Fraction of reads')





