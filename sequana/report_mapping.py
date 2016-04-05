# Import -----------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
from biokit.stats import mixture
import running_median

# Class ------------------------------------------------------------------------

class Bed_genomecov(object):
    """ Create pandas dataframe of bed file provided by bedtools genomecov (-d).
    
    :param input_filename: the input data with results of a bedtools genomecov
                           run.

    """
    def __init__(self, input_filename):
        try:
            self.df = pd.read_table(input_filename, header=None)
            self.df = self.df.rename(columns={0: "chr", 1: "pos", 2: "cov"})
        except IOError as e:
            print("I/0 error({0}): {1}".format(e.errno, e.strerror))

    def __str__(self):
        return self.df.__str__()

    def moving_average(self, n, label="ma"):
        """ Do moving average of reads coverage and create a column called 'ma'
        in data frame with results.

        :param n: window's size.

        """
        ret = np.cumsum(np.array(self.df["cov"]), dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ma = ret[n - 1:] / n
        mid = int(n / 2)
        self.df["ma"] = pd.Series(ma, index=np.arange(start=mid,
            stop=(len(ma) + mid)))

    def running_median(self, n, label="rm"):
        """ Do running median of reads coverage and create a column called 'rm'
        in data frame withe results.

        :param n: window's size.

        """
        rm = list(running_median.RunningMedian(n, self.df["cov"]))
        mid = int(n / 2)
        self.df[label] = pd.Series(rm, index=np.arange(start=mid,
            stop=(len(rm) + mid)))

    def coverage_scaling(self, method="rm", label="scale"):
        """ Normalize data with moving average of coverage and create a column 
        called 'scale' in data frame with results.
        Needs result of moving_average().

        """
        try:
            self.df[label] =  self.df["cov"] / self.df[method] 
        except KeyError:
            print("Column " + method + "is missing.\n"
                   "You must run moving_average() function before this.\n\n"
                   "Usage:\n"
                   "> mydata = Bed_genomecov('exemple.txt')\n"
                   "> mydata.moving_average(n=1000)\n"
                   "> mydata.coverage_scaling()")
            return
        self.df = self.df.replace(np.inf , np.nan)

    def _get_best_gaussian(self, results):
        diff = 100
        for i, value in enumerate(results.mus):
            if(abs(value - 1) < diff):
                diff = value
                indice = i
        return indice

    def compute_zscore(self, k=2, label="zscore", step=10):
        """ Compute zscore of coverage. 
        Needs result of coverage_scaling().

        :param k: Number gaussian predicted in mixture (default = 2)
        :param method: Column name with result of running median
        :param step: (default = 10)

        """
        try:
            mf = mixture.GaussianMixtureFitting(
                    self.df["scale"].dropna()[::step],k=k)
        except KeyError:
            print("Column '", method,"' is missing in data frame.\n"
                  "You must run coverage_scaling() function before this.\n\n"
                  "Usage:\n"
                  "> mydata = Bed_genomecov('exemple.txt')\n"
                  "> mydata.running_median(n=1000, label='rm')\n"
                  "> mydata.coverage_scaling(method='rm')\n"
                  "> mydata.compute_zscore()")
            return
        mf.estimate()
        self.gaussian = mf.results
        i = self._get_best_gaussian(mf.results)
        self.df[label] = (self.df["scale"] - mf.results["mus"][i]) / \
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

if __name__ == "__main__":
    mydata = Bed_genomecov("~/Documents/pasteur/py_dev/mapping_stats/output.txt")
    mydata.moving_average(n=30001)
    mydata.coverage_scaling()
    mydata.compute_zscore(k=2)
    plot(mydata.df["pos"], mydata.df["cov"], label="coverage")
    mydata.moving_average(n=1001)
    plot(mydata.df["pos"], mydata.df["ma"], label="w1001")
    mydata.moving_average(n=2001)
    plot(mydata.df["pos"], mydata.df["ma"], label="w2001")
    mydata.moving_average(n=5001)
    plot(mydata.df["pos"], mydata.df["ma"], label="w5001")
    mydata.moving_average(n=10001)
    plot(mydata.df["pos"], mydata.df["ma"], label="w10001")
    mydata.moving_average(n=20001)
    plot(mydata.df["pos"], mydata.df["ma"], label="w20001")
    mydata.moving_average(n=30001)
    plot(mydata.df["pos"], mydata.df["ma"], label="w30001")
    legend()

    mydata.running_median(n=1001, label="rm_1001")
    mydata.running_median(n=2001, label="rm_2001")
    mydata.running_median(n=5001, label="rm_5001")
    mydata.running_median(n=10001, label="rm_10001")
    mydata.running_median(n=30001, label="rm_30001")
    mydata.running_median(n=100001, label="rm_100001")
    plot(mydata.df["pos"], mydata.df["cov"], label="coverage")
    plot(mydata.df["pos"], mydata.df["rm_1001"], label="w1001")
    plot(mydata.df["pos"], mydata.df["rm_2001"], label="w2001")
    plot(mydata.df["pos"], mydata.df["rm_5001"], label="w5001")
    plot(mydata.df["pos"], mydata.df["rm_10001"], label="w10001")
    plot(mydata.df["pos"], mydata.df["rm_30001"], label="w30001")
    plot(mydata.df["pos"], mydata.df["rm_100001"], label="w100001")
    legend()
