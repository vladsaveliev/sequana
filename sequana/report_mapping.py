# Import -----------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
from biokit.stats import mixture

# Class ------------------------------------------------------------------------

class Bed_genomecov(object):
    """ Create pandas dataframe of bed file provided by bedtools genomecov (-d).
    
    :param input_filename: the input data with results of a bedtools genomecov
                           run.

    """
    def __init__(self, input_filename):
        try:
            self.df = pd.read_table(input_filename, header=None)
        except IOError as e:
            print("I/0 error({0}): {1}".format(e.errno, e.strerror))

    def moving_average(self, n):
        """ Do moving average of reads coverage and create a column called 'ma'
        in data frame with results.

        :param n: window's size.

        """
        ret = np.cumsum(np.array(self.df[2]), dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ma = ret[n - 1:] / n
        mid = int(n / 2)
        self.df["ma"] = pd.Series(ma, index=np.arange(start=mid,
            stop=(len(ma) + mid)))

    def coverage_scalling(self):
        """ Normalize data with moving average of coverage.
        Needs result of moving_average().

        """
        self.df["scale"] =  self.df[2] / self.df["ma"]

    def _get_best_gaussian(self, results):
        diff = 100
        for i, value in enumerate(results.mus):
            if(abs(value - 1) < diff):
                diff = value
                indice = i
        return indice

    def compute_zscore(self, k=2):
        """ Compute zscore of coverage. 
        Needs result of coverage_scaling().

        :param k: Number gaussian predicted in mixture (default = 2)

        """
        mf = mixture.GaussianMixtureFitting(self.df["scale"].dropna(), k=k)
        mf.estimate()
        self.gaussian = mf.results
        i = self._get_best_gaussian(mf.results)
        self.df["zscore"] = (ma_result - mf.results["mus"][i]) / \
            mf.results["sigmas"][i]


    def get_low_coverage(self, threshold=-3):
        try:
            return self.df.loc[self.df[3] < threshold]
        except KeyError:
            print(""" You must run compute_zscore before get low coverage.
                    self.compute_zscore(k=2, size=500)
                    """)

    def get_high_coverage(self, threshold=3):
        try:
            return self.df.loc[self.df[3] > threshold]
        except KeyError:
           print(""" You must run compute_zscore before get high coverage.
                    self.compute_zscore(k=2, size=500)
                    """)



