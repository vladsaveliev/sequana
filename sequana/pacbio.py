import pysam


class PacBioInputBAM(object):
    """PacBio utilities

    Downsample PacBio base-call BAM file 

    """
    def __init__(self, filename):
        self.filename = filename
        self.data = pysam.AlignmentFile(filename, check_sq=False)

    def reset(self):
        self.data.close()
        self.data = pysam.AlignmentFile(self.filename, check_sq=False)

    def stride(self, output_filename, stride=10):
        with pysam.AlignmentFile(output_filename,  "wb", template=self.data) as fh:

            for i, read in enumerate(self.data):
                if i % stride == 0: 
                    fh.write(read)

