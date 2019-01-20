from sequana import BAM
import glob
import pandas as pd
import pylab


class LAA():
    def __init__(self, where="bc*"):
        self.filenames = glob.glob(where + "/" + "amplicon_*summary.csv")
        self.data = [pd.read_csv(this) for this in self.filenames]

    def hist_amplicon(self, fontsize=12):
        data = [len(x) for x in self.data]
        pylab.hist(data, bins=max(data), ec="k")
        pylab.ylabel("#", fontsize=fontsize)
        pylab.ylabel("Number of amplicons per barcode", fontsize=fontsize)


class LAA_Assembly():
    """

    Input is a SAM/BAM from the mapping of amplicon onto a known reference.
    Based on the position, we can construct the new reference.

    """
    def __init__(self, filename):
        self.bam = BAM(filename)


    def build_reference(self):
        self.bam.reset()
        # scan BAM file assuming it is small
        aa = [a for a in self.bam]

        # retrieve data of interest
        data = [(a.pos, {
                    "name":a.query_name,
                    "sequence": a.query_sequence,
                    "cigar": a.cigarstring,
                    "position": a.pos,
                    "qstart": a.qstart,
                    "qend": a.qend}) for a in aa]

        # sort by starting position
        data.sort(key=lambda x: x[0])

        for i, read in enumerate(data):
            read = read[1]
            if i == 0:
                sequence = read["sequence"]     # 2 is query_sequence
            else:
                pr = data[i-1][1]   # previous read
                L = len(pr["sequence"])
                end_position_pr = pr['position'] - pr['qstart'] + L 

                # overlap between previous read and this one
                overlap = end_position_pr - (read['position'] - read['qstart']) +0
                print(overlap)
                print(pr['position'], pr['qstart'], L, end_position_pr)
                print(read['position'], read['qstart'])
                sequence = sequence + read["sequence"][overlap+1:]

        # argmax([sum(a==b for a,b in zip(X[-i:] , Y[:i]))/float(i+1) for i in range(1000)])
        return sequence

    def save_fasta(self, filename, sequence=None):
        if sequence is None:
            sequence = self.build_reference()


        with open(filename, "w") as fout:
            fout.write(">test\n{}".format(sequence))


