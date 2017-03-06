
from sequana.lazy import pylab
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
import collections #lazy ?
import pysam
from biokit.viz import hist2d

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
        self._df = None
        self._nb_pass = None

    def __len__(self):
        if self._N is None:
            df = self._get_df()
        return self._N

    def __str__(self):
        return "Length: {}".format(len(self))

    def _get_df(self):
        if self._df is None:
            self.reset()
            N = 0
            
            all_results = []
            for read in self.data:
                res = []
                # count reads
                N += 1
                if (N % 10000) == 0:
                    print("Read %d sequences" %N)
                #res[0] = read length
                res.append(read.query_length)
                # res[1] = GC content
                c = collections.Counter(read.query_sequence)
                res.append( (c['g'] + c['G'] + c['c'] + c['C'])/float(sum(c.values())) )
                # res[2] = snr A
                # res[3] = snr C
                # res[4] = snr G
                # res[5] = snr T
                snr = list([x for x in read.tags if x[0]=='sn'][0][1])
                res = res + snr
                #res[6] = ZMW name
                res.append(read.qname.split('/')[1])
                
                # aggregate results
                all_results.append(res)

            self._df = pd.DataFrame(all_results, columns=['read_length','GC_content','snr_A','snr_C','snr_G','snr_T','ZMW'])
            self._N = N
            self.reset()     
        return self._df

    df = property(_get_df)

    def _get_ZMW_passes(self):
        print()
        if self._nb_pass is None:
            if self._df is None:
                self._get_df()

            zmw_passes = collections.Counter(self._df.loc[:,'ZMW'])
            distrib_nb_passes = [zmw_passes[z] for z in zmw_passes.keys()]
            self._nb_pass = collections.Counter(distrib_nb_passes)
        return self._nb_pass

    nb_pass = property(_get_ZMW_passes)

    def reset(self):
        self.data.close()
        self.data = pysam.AlignmentFile(self.filename, check_sq=False)

    def stride(self, output_filename, stride=10):
        self.reset()
        with pysam.AlignmentFile(output_filename,  "wb", template=self.data) as fh:

            for i, read in enumerate(self.data):
                if i % stride == 0: 
                    fh.write(read)


    def hist_snr(self, bins=50, alpha=0.5, hold=False, fontsize=12,
                grid=True,xlabel="SNR",ylabel="#"):
        """Plot histogram of the ACGT SNRs for all reads"""
        if self._df is None:
            self._get_df()

        if hold is False:
            pylab.clf()
        pylab.hist(self._df.loc[:,'snr_A'], alpha=alpha, label="A", bins=bins)
        pylab.hist(self._df.loc[:,'snr_C'], alpha=alpha, label="C", bins=bins)
        pylab.hist(self._df.loc[:,'snr_G'], alpha=alpha, label="G", bins=bins)
        pylab.hist(self._df.loc[:,'snr_T'], alpha=alpha, label="T", bins=bins)
        pylab.legend()
        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        if grid is True:
            pylab.grid(True)

    def hist_ZMW_subreads(self, hold=False, fontsize=12,
                            grid=True,xlabel="Number of ZMW passes",ylabel="#"):
        """
        Plot histogram of number of reads per ZMW
        """
        if self._nb_pass is None:
            self._get_ZMW_passes()

        max_nb_pass = max(self._nb_pass.keys())
        k = range(1,max_nb_pass+1)
        val = [self._nb_pass[i] for i in k]

        # histogram nb passes
        if hold is False:
            pylab.clf()
        pylab.hist(k, weights=val, bins=max_nb_pass)
        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        pylab.yscale('log')
        pylab.title("Number of ZMW passes",fontsize=fontsize)
        if grid is True:
            pylab.grid(True)

    def hist_GC(self, bins=50, hold=False, fontsize=12,
                grid=True,xlabel="GC %",ylabel="#"):
        """Plot histogram GC content"""

        if self._df is None:
            self._get_df()
        mean_GC =  np.mean(self._df.loc[:,'GC_content'])

        # histogram GC percent
        if hold is False:
            pylab.clf()
        pylab.hist(self._df.loc[:,'GC_content'], bins=bins)
        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        pylab.title("GC %%  \n Mean GC : %.2f" %(mean_GC), fontsize=fontsize)
        if grid is True:
            pylab.grid(True)

    def plot_GC_read_len(self, alpha=0.07, hold=False, fontsize=12,
                grid=True,xlabel="GC %",ylabel="#"):
        """Plot GC content versus read length"""

        if self._df is None:
            self._get_df()
        mean_len =  np.mean(self._df.loc[:,'read_length'])
        mean_GC =  np.mean(self._df.loc[:,'GC_content'])

        if hold is False:
            pylab.clf()
        data = self._df.loc[:,['read_length','GC_content']]
        h = hist2d.Hist2D(data)
        res = h.plot(bins=[40,40], contour=False, nnorm='log', Nlevels=6)
        #pylab.plot(self._df.loc[:,'read_length'] , self._df.loc[:,'GC_content'], 'bo', alpha=alpha)
        pylab.xlabel("Read length", fontsize=12)
        pylab.ylabel("GC %", fontsize=12)
        pylab.title("GC % vs length \n Mean length : %.2f , Mean GC : %.2f" %(mean_len, mean_GC))





