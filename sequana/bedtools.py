# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Utilities for the genome coverage"""
import os

import pandas as pd
import numpy as np
import pylab
from biokit.stats import mixture

from sequana import running_median


class Genomecov(object):
    """Create a dataframe of BED file provided by bedtools genomecov (-d)


    Example:

    .. plot::
        :include-source:

        from sequana import Genomecov, sequana_data
        filename = sequana_data("test_bedcov.bed", "testing")

        gencov = Genomecov(filename)
        gencov.running_median(n=1001 )
        gencov.coverage_scaling()
        gencov.compute_zscore()
        gencov.plot_coverage()

    Results are stored in a dataframe named :attr:`df`.

    """
    def __init__(self, input_filename=None):
        """.. rubric:: constructor

        :param str input_filename: the input data with results of a bedtools
            genomecov run. This is just a 3-column file. The first column is a
            string, second column is the base postion and third is the coverage.

        """
        try:
            self.df = pd.read_table(input_filename, header=None)
            self.df = self.df.rename(columns={0: "chr", 1: "pos", 2: "cov"})
        except IOError as e:
            print("I/0 error({0}): {1}".format(e.errno, e.strerror))

    def __str__(self):
        return self.df.__str__()

    def __len__(self):
        return self.df.__len__()

    def moving_average(self, n, label="ma"):
        """Compute moving average of reads coverage

        :param n: window's size.

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *ma*.

        """
        ret = np.cumsum(np.array(self.df["cov"]), dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ma = ret[n - 1:] / n
        mid = int(n / 2)
        self.df["ma"] = pd.Series(ma, index=np.arange(start=mid,
            stop=(len(ma) + mid)))

    def running_median(self, n, circular=False):
        """Compute running median of reads coverage

        :param int n: window's size.
        :param bool circular: if a mapping is circular (e.g. bacteria
            whole genome sequencing), set to True

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *rm*.

        """
        mid = int(n / 2)# in py2/py3 the division (integer or not) has no impact
        if circular:
            cov = list(self.df["cov"])
            cov = cov[-mid:] + cov + cov[:mid]
            rm = running_median.RunningMedian(cov, n).run()
            self.df["rm"] = pd.Series(rm)
        else:
            rm = running_median.RunningMedian(self.df["cov"], n).run()
            self.df["rm"] = pd.Series(rm)
            #, index=np.arange(start=mid,
            #    stop=(len(rm) + mid)))

    def coverage_scaling(self):
        """Normalize data with moving average of coverage

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *scale*.

        .. note:: Needs to call :meth:`running_median`
        """
        try:
            self.df["scale"] =  self.df["cov"] / self.df["rm"]
        except KeyError:
            print("Column rm (running median) is missing.\n",
                    self.__doc__)
            return
        self.df = self.df.replace(np.inf , np.nan)

    def _get_best_gaussian(self, results):
        diff = 100
        for i, value in enumerate(results.mus):
            if abs(value - 1) < diff:
                diff = abs(value - 1)
                indice = i
        return {"mu": results.mus[indice], "sigma": results.sigmas[indice]}

    def compute_zscore(self, k=2, step=10):
        """ Compute zscore of coverage.

        :param k: Number gaussian predicted in mixture (default = 2)
        :param step: (default = 10)

        .. note:: needs to call :meth:`coverage_scaling` before hand.

        """
        try:
            self.mixture_fitting = mixture.GaussianMixtureFitting(
                    self.df["scale"].dropna()[::step],k=k)
        except KeyError:
            print("Column 'scale' is missing in data frame.\n"
                  "You must scale your data with coverage_scaling\n\n",
                  self.__doc__)
            return
        self.mixture_fitting.estimate()
        self.gaussians = self.mixture_fitting.results
        self.best_gaussian = self._get_best_gaussian(self.mixture_fitting.results)
        self.df["zscore"] = (self.df["scale"] - self.best_gaussian["mu"]) / \
            self.best_gaussian["sigma"]

    def get_low_coverage(self, threshold=-3, start=None, stop=None):
        """Keep position with zscore lower than INT and return a data frame.

        :param int threshold: on the zscore
        :param int start: lower bound to select a subset of the data
        :param int stop:  upper bound to select a subset of the data
        :return: a dataframe from :class:`FilteredGenomecov`
        """
        try:
            return FilteredGenomecov(self.df[start:stop].loc[self.df["zscore"]
                < threshold])
        except KeyError:
            print("Column zscore is missing in data frame.\n"
                  "You must run compute_zscore before get low coverage.\n\n",
                  self.__doc__)

    def get_high_coverage(self, threshold=3, start=None, stop=None):
        """Keep position with zscore higher than INT and return a data frame.

        :param int threshold: on the zscore
        :param int start: lower bound to select a subset of the data
        :param int stop:  upper bound to select a subset of the data
        :return: a dataframe from :class:`FilteredGenomecov`
        """
        try:
            return FilteredGenomecov(self.df[start:stop].loc[self.df["zscore"]
                > threshold])
        except KeyError:
            print("Column zscore is missing in data frame.\n"
                  "You must run compute_zscore before get low coverage.\n\n",
                    self.__doc__)

    def plot_coverage(self, fontsize=16, filename=None,
            low_threshold=-3, high_threshold=3):
        """ Plot coverage as a function of base position.

        In addition, the running median and coverage confidence corrsponding to
        the lower and upper  zscore thresholds

        """
        high_zcov = (high_threshold * self.best_gaussian["sigma"] +
                self.best_gaussian["mu"]) * self.df["rm"]
        low_zcov = (low_threshold * self.best_gaussian["sigma"] +
                self.best_gaussian["mu"]) * self.df["rm"]

        pylab.clf()
        pylab.xlim(0,self.df["pos"].iloc[-1])
        p1, = pylab.plot(self.df["cov"], color="b", label="Coverage",
                linewidth=1)
        p2, = pylab.plot(self.df["rm"], color="r", linewidth=1,
                label="Running median")
        p3, = pylab.plot(high_zcov, linewidth=1, color="r", ls="--",
                label="Thresholds")
        p4, = pylab.plot(low_zcov, linewidth=1, color="r", ls="--")

        pylab.legend([p1, p2, p3], [p1.get_label(), p2.get_label(),
                p3.get_label()], loc="best")
        pylab.xlabel("Position", fontsize=fontsize)
        pylab.ylabel("Coverage", fontsize=fontsize)
        try:
            pylab.tight_layout()
        except:
            pass
        if filename:
            pylab.savefig(filename)

    def plot_hist(self, fontsize=16, bins=100,filename=None ):
        """ Barplot of zscore

        """
        pylab.clf()
        self.df["zscore"].hist(grid=True, color="b", bins=100)
        pylab.xlabel("Z-Score", fontsize=fontsize)
        try:
            pylab.tight_layout()
        except:
            pass
        if filename:
            pylab.savefig(filename)

    def write_csv(self, filename, start=None, stop=None, header=True):
        """ Write CSV file of the dataframe.

        :param filename: csv output filename.
        :param header: boolean which determinate if the header is written.

        """
        try:
            labels=["pos", "cov", "rm"]
            self.df[labels][start:stop].to_csv(filename, header=header)
        except NameError:
            print("You must set the file name")
        except KeyError:
            print("Labels doesn't exist in the data frame")


class FilteredGenomecov(object):
    """Class used within :class:`Genomecov` to select a subset of the original Genomecov

    :target: developers only
    """
    def __init__(self, df):
        """ .. rubric:: constructor

        :param df: dataframe with filtered position used within
            :class:`Genomecov`. Must contain the following columns:
            ["pos", "cov", "rm", "zscore"]

        """
        self.df = df

    def __str__(self):
        return self.df.__str__()

    def __len__(self):
        return self.df.__len__()

    def _merge_row(self, start, stop):
        chrom = self.df["chr"][start]
        cov = np.mean(self.df["cov"].loc[start:stop])
        rm = np.mean(self.df["rm"].loc[start:stop])
        zscore = np.mean(self.df["zscore"].loc[start:stop])
        size = stop - start
        return {"chr": chrom, "start": start + 1, "stop": stop, "size": size,
                "mean_cov": cov, "mean_rm": rm, "mean_zscore": zscore}

    def merge_region(self, threshold, zscore_label="zscore"):
        """Merge position side by side of a data frame.

        .. todo:: to be documented
        """
        flag = False
        start = 1
        stop = 1
        prev = 1

        merge_df = pd.DataFrame(columns=["chr", "start", "stop", "size",
            "mean_cov", "mean_rm", "mean_zscore"])
        int_column = ["start", "stop", "size"]
        merge_df[int_column] = merge_df[int_column].astype(int)

        for pos, zscore in zip(self.df["pos"], self.df[zscore_label]):
            stop = pos
            if stop - 1 == prev:
                prev = stop
            else:
                if flag:
                    merge_df = merge_df.append(self._merge_row(start - 1, prev),
                        ignore_index=True)
                    flag = False
                start = stop
                prev = stop
            if abs(zscore) > abs(threshold):
                flag = True

        if start < stop and flag:
            merge_df = merge_df.append(self._merge_row(start - 1, prev),
                    ignore_index=True)
        return merge_df
