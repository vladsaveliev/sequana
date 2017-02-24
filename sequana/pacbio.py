
from sequana.lazy import pylab
import collections #lazy ?
from sequana.lazy import numpy as np
import pysam


class PacBioInputBAM(object):
    """PacBio utilities

    Downsample PacBio base-call BAM file 

    TODO:

        number of sub reads per ZMW > hist_ZMW_subreads(self)

    """
    def __init__(self, filename):
        self.filename = filename
        self.data = pysam.AlignmentFile(filename, check_sq=False)
        self._N = None
        self.GC_percent = None
        self.all_read_len = None
        self.nb_pass = None

    def __len__(self):
        if self._N is None:
            self.reset()
            self._N = sum([1 for this in self.data])
            self.reset()
        return self._N

    def __str__(self):
        return "Length: {}".format(len(self))

    def reset(self):
        self.data.close()
        self.data = pysam.AlignmentFile(self.filename, check_sq=False)

    def stride(self, output_filename, stride=10):
        self.reset()
        with pysam.AlignmentFile(output_filename,  "wb", template=self.data) as fh:

            for i, read in enumerate(self.data):
                if i % stride == 0: 
                    fh.write(read)

    def GC_content(self):
        """GC content for each read
        Save result of GC content in self.GC_percent if not already present"""
        if self.GC_percent is None:
            self.reset()
            GC_percent = []
            for read in self.data:
                c = collections.Counter(read.query_sequence)
                GC_percent.append( (c['g'] + c['G'] + c['c'] + c['C'])/float(sum(c.values())) )
            self.reset()
            self.GC_percent = GC_percent

    def read_len(self):
        """Read len for each read)"""
        if self.all_read_len is None:
            self.reset()
            self.all_read_len = [read.query_length for read in self.data]

    def hist_snr(self, bins=50, alpha=0.5, hold=False, fontsize=12,
                grid=True,xlabel="SNR",ylabel="#"):
        """Plot histogram of the ACGT SNRs for all reads"""
        self.reset()
        reads = [[x for x in read.tags if x[0]=='sn'] for read in self.data]
        snr = [x[0][1] for x in reads]
        if hold is False:
            pylab.clf()
        pylab.hist([x[0] for x in snr], alpha=alpha, label="A", bins=bins)
        pylab.hist([x[1] for x in snr], alpha=alpha, label="C", bins=bins)
        pylab.hist([x[2] for x in snr], alpha=alpha, label="G", bins=bins)
        pylab.hist([x[3] for x in snr], alpha=alpha, label="T", bins=bins)
        pylab.legend()
        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        if grid is True:
            pylab.grid(True)

    def hist_ZMW_subreads(self):
        """
        Plot histogram of number of reads per ZMW
        Save counter in nb_pass
        """
        if self.nb_pass is None:
            self.reset()
            zmw = [read.qname.split('/')[1] for read in self.data]
            zmw_passes = collections.Counter(zmw)
            distrib_nb_passes = [zmw_passes[z] for z in zmw_passes.keys()]
            max_nb_pass = max(distrib_nb_passes)
            self.nb_pass = collections.Counter(distrib_nb_passes)

        max_nb_pass = max(self.nb_pass.keys())
        k = range(1,max_nb_pass+1)
        print(k)
        val = [self.nb_pass[i] for i in k]
        print(val)
        # histogram nb passes
        pylab.hist(k, weights=val, bins=max_nb_pass)
        pylab.xlabel("Nb passes")
        pylab.ylabel("ZMW")
        pylab.yscale('log')
        pylab.title("Nb passes")

    def hist_GC(self, bins=50):
        """Plot histogram GC content"""
        self.GC_content()
        mean_GC =  np.mean(self.GC_percent)

        # histogram GC percent
        pylab.hist(self.GC_percent, bins=bins)
        pylab.xlabel("GC percent")
        pylab.title("GC content  \n Mean GC : %.2f" %(mean_GC))

    def plot_GC_read_len(self, alpha=0.07):
        """Plot GC content versus read length"""
        self.read_len()
        self.GC_content()
        mean_len =  np.mean(self.all_read_len)
        mean_GC =  np.mean(self.GC_percent)

        pylab.plot(self.all_read_len , self.GC_percent, 'bo', alpha=alpha)
        pylab.xlabel("Read length")
        pylab.ylabel("GC percent")
        pylab.title("GC content vs length \n Mean length = %.2f , Mean GC : %.2f" %(mean_len, mean_GC))





