"""Utilities to manipulate FASTQ and Reads

"""
import io
import time
import zlib
from itertools import islice
import gzip
import subprocess
from functools import wraps

import numpy as np
import pandas as pd
from easydev import do_profile, Progress
import pylab

import pysam
try:
    from itertools import izip_longest
except:
    from itertools import zip_longest as izip_longest

from sequana.fastq import FastQ
from pysam import FastxFile

# for filter fastq files. see below in FastQ for the usage
# we want to take 4 lines at a time (assuming there is no empty lines)

__all__ = ["FastA"]

# cannot inherit from FastxFile (no object in the API ?)
class FastA(object):
    """Class to handle FastA files


    """
    def __init__(self, filename, verbose=False):
        self.fasta = FastxFile(filename)
        self.N = len([x for x in FastxFile(filename)])

    def __iter__(self):
        return self

    def __next__(self): # python 3
        return self.next()

    def next(self): # python 2
        # reads 4 lines
        try:
            d = next(self.fasta)
            return d
        except KeyboardInterrupt:
            # THis should allow developers to break a loop that takes too long
            # through the read to run forever
            self.fasta.close()
            self.fasta = FastxFile(self.fasta.filename)
        except:
            self.fasta.close()
            self.fasta = FastxFile(self.fasta.filename)
            raise StopIteration

        return d

    def to_fasta(self):
        pass

    def __len__(self):
        return self.N
