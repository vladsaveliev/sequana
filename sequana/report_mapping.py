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

    def __str__(self):
        return self.df.__str__()

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

    def coverage_scaling(self):
        """ Normalize data with moving average of coverage and create a column 
        called 'scale' in data frame with results.
        Needs result of moving_average().

        """
        try:
            self.df["scale"] =  self.df[2] / self.df["ma"]
        except KeyError:
            print("Column 'ma' is missing.\n"
                   "You must run moving_average() function before this.\n\n"
                   "Usage:\n"
                   "> mydata = Bed_genomecov('exemple.txt')\n"
                   "> mydata.moving_average(n=1000)\n"
                   "> mydata.coverage_scaling()")
            return



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
        try:
            mf = mixture.GaussianMixtureFitting(self.df["scale"].dropna(), k=k)
        except KeyError:
            print("Column 'scale' is missing in data frame.\n"
                  "You must run coverage_scaling() function before this.\n\n"
                  "Usage:\n"
                  "> mydata = Bed_genomecov('exemple.txt')\n"
                  "> mydata.moving_average(n=1000)\n"
                  "> mydata.coverage_scaling()\n"
                  "> mydata.compute_zscore()")
            return
        mf.estimate()
        self.gaussian = mf.results
        i = self._get_best_gaussian(mf.results)
        self.df["zscore"] = (self.df["scale"] - mf.results["mus"][i]) / \
            mf.results["sigmas"][i]


    def get_low_coverage(self, threshold=-3):
        """Keep position with zscore lower than INT and return a data frame.

        :param threshold: Integer
        """
        try:
            return self.df.loc[self.df["zscore"] < threshold]
        except KeyError:
            print("Column 'zscore' is missing in data frame.\n"
                  "You must run compute_zscore before get low coverage.\n\n"
                  "Usage:\n"
                  "> mydata = Bed_genomecov('exemple.txt')\n"
                  "> mydata.moving_average(n=1000)\n"
                  "> mydata.coverage_scaling()\n"
                  "> mydata.compute_zscore(k=2)")

    def get_high_coverage(self, threshold=3):
        """Keep position with zscore higher than INT and return a data frame.

        :param threshold: Integer
        """
        try:
            return self.df.loc[self.df["zscore"] > threshold]
        except KeyError:
            print("Column 'zscore' is missing in data frame.\n"
                  "You must run compute_zscore before get low coverage.\n\n"
                  "Usage:\n"
                  "> mydata = Bed_genomecov('exemple.txt')\n"
                  "> mydata.moving_average(n=1000)\n"
                  "> mydata.coverage_scaling()\n"
                  "> mydata.compute_zscore(k=2)")



