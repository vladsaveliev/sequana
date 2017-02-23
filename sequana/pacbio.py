
from sequana.lazy import pylab
import pysam


class PacBioInputBAM(object):
    """PacBio utilities

    Downsample PacBio base-call BAM file 

    TODO:

        number of sub reads per ZMW

    SNR ACGT:

        reads = [[x for x in read.tags if x[0]=='sn'] for read in b.data]
        snr = [x[0][1] for x in reads]
        hist([x[0] for x in snr], alpha=0.5, label="A", bins=50)
        hist([x[1] for x in snr], alpha=0.5, label="C", bins=50)
        hist([x[2] for x in snr], alpha=0.5, label="G", bins=50)
        hist([x[3] for x in snr], alpha=0.5, label="T", bins=50)
        legend()

    """
    def __init__(self, filename):
        self.filename = filename
        self.data = pysam.AlignmentFile(filename, check_sq=False)
        self._N = None

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

    def hist_snr(self, bins=50, alpha=0.5):
        """Plot histogram of the ACGT SNRs for all reads"""
        self.reset()
        reads = [[x for x in read.tags if x[0]=='sn'] for read in self.data]
        snr = [x[0][1] for x in reads]
        pylab.hist([x[0] for x in snr], alpha=alpha, label="A", bins=bins)
        pylab.hist([x[1] for x in snr], alpha=alpha, label="C", bins=bins)
        pylab.hist([x[2] for x in snr], alpha=alpha, label="G", bins=bins)
        pylab.hist([x[3] for x in snr], alpha=alpha, label="T", bins=bins)
        pylab.legend()



