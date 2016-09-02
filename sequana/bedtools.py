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
import sys

import pandas as pd
import numpy as np
import pylab
from biokit.stats import mixture

from sequana import running_median
from sequana.tools import gc_content
from sequana.tools import genbank_features_parser
from easydev import TempFile

class GenomeCov(object):
    """Create a dataframe list of BED file provided by bedtools genomecov (-d)


    Example:

    .. plot::
        :include-source:

        from sequana import GenomeCov, sequana_data
        filename = sequana_data("test_bedcov.bed", "testing")

        gencov = GenomeCov(filename)
        for chrom in gencov:
            chrom.running_median(n=3001)
            chrom.compute_zscore()
            chrom.plot_coverage()
        gencov[0].plot_coverage()

    Results are stored in a list of :class:`ChromosomeCov` named :attr:`chr_list`.

    """

    def __init__(self, input_filename=None):
        """.. rubric:: constructor

        :param str input_filename: the input data with results of a bedtools
            genomecov run. This is just a 3-column file. The first column is a
            string, second column is the base postion and third is the coverage.

        """

        df = pd.read_table(input_filename, header=None)
        try:
            df = df.rename(columns={0: "chr", 1: "pos", 2: "cov"})
            df = df.set_index("chr", drop=False)
            # Create a list of ChromosomeCov for each chromosome present in the
            # bedtools.
            self.chr_list = [ChromosomeCov(df.loc[key]) for key in
                    df.index.unique()]
            ChromosomeCov.count = 0
        except IOError as e:
            print("I/0 error({0}): {1}".format(e.errno, e.strerror))

        # Set the link to this instance in each chromosome
        # useful if one wants to recompute GC content with different window
        for chrom in self.chr_list:
            chrom.bed = self

    def __getitem__(self, index):
        return self.chr_list[index]

    def __iter__(self):
        return self.chr_list.__iter__()

    def compute_gc_content(self, fasta_file, window_size=101, circular=False):
        """ Compute GC content of genome sequence.

        :param str fasta_file: fasta file name.
        :param int window_size: size of the sliding window.
        :param bool circular: if the genome is circular (like bacteria
            chromosome)

        Store the results in the :attr:`ChromosomeCov.df` attribute (dataframe)
            with a column named *gc*.

        """
        gc_dict = gc_content(fasta_file, window_size, circular)
        for chrom in self.chr_list:
            chrom.df["gc"] = gc_dict[chrom.chrom_name]
            chrom._ws_gc = window_size

    def get_stats(self):
        stats = {}
        for chrom in self.chr_list:
            stats[chrom.chrom_name] = chrom.get_stats()
        return stats


class ChromosomeCov(object):
    """Class used within :class:`GenomeCov` to select a chromosome of the
    original GenomeCov.

    Example:

    .. plot::
        :include-source:

        from sequana import GenomeCov, sequana_data
        filename = sequana_data("test_bedcov.bed", "testing")

        gencov = GenomeCov(filename)

        chrcov = gencov[0]
        chrcov.running_median(n=3001)
        chrcov.compute_zscore()
        chrcov.plot_coverage()

    Results are stored in a dataframe named :attr:`df`.

    """
    count = 0

    def __init__(self, df=None):
        """.. rubric:: constructor

        :param df: dataframe with position for a chromosome used within
            :class:`GenomeCov`. Must contain the following columns:
            ["chr", "pos", "cov"]

        """
        # Chromosome position becomes the index
        ChromosomeCov.count += 1
        self.chrom_index = ChromosomeCov.count
        self.df = df.set_index("pos", drop=False)
        self.chrom_name = df["chr"].iloc[0]

    def __str__(self):
        stats = self.get_stats()
        BOC = stats['BOC']
        CV = stats['CV']
        txt = "\nGenome length: %s" % int(len(self.df))
        txt += "\nSequencing depth (DOC): %8.2f " % stats['DOC']
        txt += "\nSequencing depth (median): %8.2f " % stats['median']
        txt += "\nBreadth of coverage (BOC): %.2f " % BOC
        txt += "\nGenome coverage standard deviation : %8.2f " % stats['std']
        txt += "\nGenome coverage coefficient variation : %8.2f " % CV
        return txt

    def __len__(self):
        return self.df.__len__()

    def get_size(self):
        return self.__len__()

    def get_mean_cov(self):
        return self.df["cov"].mean()

    def get_var_coef(self):
        return np.sqrt(self.df["cov"].var()) / self.get_mean_cov()

    def moving_average(self, n, circular=False):
        """Compute moving average of reads coverage

        :param n: window's size.

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *ma*.

        """
        N = len(self.df['cov'])
        assert n < N/2
        from sequana.stats import moving_average

        ret = np.cumsum(np.array(self.df["cov"]), dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        ma = ret[n - 1:] / n
        mid = int(n / 2)
        self.df["ma"] = pd.Series(ma, index=np.arange(start=mid,
            stop=(len(ma) + mid)))

        if circular:
            # FIXME: shift of +-1 as compared to non circular case...
            # shift the data and compute the moving average
            self.data = list(self.df['cov'].values[N-n:]) +\
                list(self.df['cov'].values) + \
                list(self.df['cov'].values[0:n])
            ma = moving_average(self.data, n)
            self.ma = ma[n//2+1:-n//2]
            self.df["ma"] = pd.Series(self.ma, index=self.df['cov'].index)

    def running_median(self, n, circular=False):
        """Compute running median of reads coverage

        :param int n: window's size.
        :param bool circular: if a mapping is circular (e.g. bacteria
            whole genome sequencing), set to True

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *rm*.

        """
        mid = int(n / 2)# in py2/py3 the division (integer or not) has no impact
        self.range = [None, None]
        cover = list(self.df["cov"])
        try:
            if circular:
                cover = cover[-mid:] + cover + cover[:mid]
                rm = running_median.RunningMedian(cover, n).run()
                self.df["rm"] = rm[mid:-mid]
            else:
                rm = running_median.RunningMedian(cover, n).run()
                self.df["rm"] = rm
                # set up slice for gaussian prediction
                self.range = [mid, -mid]
        except:
            self.df["rm"] = self.df["cov"]

    def get_evenness(self):
        """Return Evenness of the coverage

        :Reference: Konrad Oexle, Journal of Human Genetics 2016, Evaulation
            of the evenness score in NGS.

        work before or after normalisation but lead to different results.

        """
        from sequana.stats import evenness
        return evenness(self.df['cov'])

    def get_cv(self):
        """coefficient variation

        defined as sigma / mu

        To get percentage, you must multiply by 100

        .. note:: should be used for ratio scale data (e.g., non negative only)
        """
        sigma = self.df['cov'].std()
        mu = self.df['cov'].mean()
        return sigma/mu

    def _coverage_scaling(self):
        """Normalize data with moving average of coverage

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *scale*.

        .. note:: Needs to call :meth:`running_median`


        """
        if "rm" not in self.df.columns:
            txt = "Column rm (running median) is missing.\n" +  self.__doc__
            print(txt)
            raise KeyError
        else:
            self.df["scale"] =  self.df["cov"] / self.df["rm"]
        self.df = self.df.replace(np.inf, np.nan)
        self.df = self.df.replace(-np.inf, np.nan)

    def _get_best_gaussian(self):
        results = self.mixture_fitting.results
        indice = np.argmax(results.pis)
        return {"mu": results.mus[indice], "sigma": results.sigmas[indice]}

    def compute_zscore(self, k=2, step=10, use_em=True):
        """ Compute zscore of coverage and normalized coverage.

        :param int k: Number gaussian predicted in mixture (default = 2)
        :param int step: (default = 10). This parameter is used to speed
            up computation and is ignored if the length of the coverage/sequence
            is below 100,000

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *zscore*.

        .. note:: needs to call :meth:`running_median` before hand.

        """
        # normalize coverage
        self._coverage_scaling()

        data = self.df['scale'][self.range[0]:self.range[1]]

        if len(data) < 100000:
            step = 1

        # remove nan and inf values
        data = data.replace(0, np.nan)
        data = data.dropna()

        if data.empty:
            data = np.full(len(self.df), 1, dtype=int)
            self.df['scale'] = data

        if use_em:
            self.mixture_fitting = mixture.EM(
                data[::step])
            self.mixture_fitting.estimate(k=k)
        else:
            self.mixture_fitting = mixture.GaussianMixtureFitting(
                data[::step],k=k)
            self.mixture_fitting.estimate()

        # keep gaussians informations
        self.gaussians = self.mixture_fitting.results
        self.best_gaussian = self._get_best_gaussian()

        # warning when sigma is equal to 0
        if self.best_gaussian["sigma"] == 0:
            print("WARNING: A problem related to gaussian prediction is "
                  "detected. Be careful, Sigma is equal to 0.")
            self.df["zscore"] = np.zeros(len(self.df), dtype=int)
        else:
            self.df["zscore"] = (self.df["scale"] - self.best_gaussian["mu"]) / \
                self.best_gaussian["sigma"]

    def get_centralness(self, threshold=3):
        r"""Proportion of central (normal) genome coverage

        assuming a 3 sigma normality.

        This is 1 - (number of non normal data) / (total length)

        .. note:: depends slightly on :math:`W` the running median window
        """
        l1 = len(self.get_low_coverage(-threshold))
        l2 = len(self.get_high_coverage(threshold))
        return 1 - (l1+l2) / float(len(self))

    def get_roi(self, first_thr=3, second_thr=1.5, features=None):
        """Keep position with zscore lower than INT and return a data frame.

        :param int first_thr: principal threshold on zscore
        :param int second_thr: secondary threshold on zscore
        :return: a dataframe from :class:`FilteredGenomeCov`
        """
        try:
            if features:
                return FilteredGenomeCov(self.df.loc[abs(self.df["zscore"]) > 
                    second_thr], first_thr, features[self.chrom_name])
            else:
                return FilteredGenomeCov(self.df.loc[abs(self.df["zscore"]) > 
                    second_thr], first_thr)
        except KeyError:
            print("Column zscore is missing in data frame.\n"
                  "You must run compute_zscore before get low coverage.\n\n",
                  self.__doc__) 

    def plot_coverage(self, filename=None, threshold=3, fontsize=16):
        """ Plot coverage as a function of base position.

        In addition, the running median and coverage confidence corresponding to
        the lower and upper  zscore thresholds

        """
        # z = (X/rm - \mu ) / sigma

        high_zcov = (threshold * self.best_gaussian["sigma"] +
                self.best_gaussian["mu"]) * self.df["rm"]
        low_zcov = (-threshold * self.best_gaussian["sigma"] +
                self.best_gaussian["mu"]) * self.df["rm"]

        pylab.clf()
        ax = pylab.gca()
        ax.set_axis_bgcolor('#eeeeee')
        pylab.xlim(0,self.df["pos"].iloc[-1])
        p1, = pylab.plot(self.df["cov"], color="k", label="Coverage",
                linewidth=1)
        p2, = pylab.plot(self.df["rm"], color="#0099cc", linewidth=1,
                label="Running median")
        p3, = pylab.plot(high_zcov, linewidth=1, color="r", ls="--",
                label="Thresholds")
        p4, = pylab.plot(low_zcov, linewidth=1, color="r", ls="--",
            label="_nolegend_")

        pylab.legend([p1, p2, p3], [p1.get_label(), p2.get_label(),
                p3.get_label()], loc="best")
        pylab.xlabel("Position", fontsize=fontsize)
        pylab.ylabel("Per-base coverage", fontsize=fontsize)
        pylab.grid(True)
        pylab.ylim([0, pylab.ylim()[1]])
        try:
            pylab.tight_layout()
        except:
            pass
        if filename:
            pylab.savefig(filename)

    def _set_bins(self, df, binwidth):
        try:
            bins = np.arange(min(df), max(df) + binwidth, binwidth)
        except ValueError:
            bins = 100
        return bins

    def plot_hist_zscore(self, fontsize=16, filename=None, max_z=6,
            binwidth=0.5, **hist_kargs):
        """ Barplot of zscore

        """
        pylab.clf()
        bins = self._set_bins(self.df["zscore"][self.range[0]:self.range[1]],
                binwidth)
        self.df["zscore"][self.range[0]:self.range[1]].hist(grid=True,
                bins=bins, **hist_kargs)
        pylab.xlabel("Z-Score", fontsize=fontsize)
        try:
            pylab.tight_layout()
        except:
            pass
        if filename:
            pylab.savefig(filename)

    def plot_hist_normalized_coverage(self, filename=None, binwidth=0.1,
            max_z=4):
        """ Barplot of normalized coverage with gaussian fitting

        """
        pylab.clf()
        # if there are a NaN -> can't set up binning
        data_scale = self.df["scale"][self.range[0]:self.range[1]].dropna()
        bins = self._set_bins(data_scale, binwidth)
        self.mixture_fitting.plot(bins=bins, Xmin=0, Xmax=max_z)
        pylab.grid(True)
        pylab.xlim([0,max_z])
        pylab.xlabel("Normalised per-base coverage")
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

    def plot_gc_vs_coverage(self, bins=None, Nlevels=6, fontsize=20, norm="log",
            ymin=0, ymax=100, contour=True, **kargs):

        if Nlevels is None or Nlevels==0:
            contour = False

        data = self.df[['cov','gc']].copy()
        data['gc'] *= 100
        data = data.dropna()
        if bins is None:
            bins = [100, min(int(data['gc'].max()-data['gc'].min()+1),
                max(5,self._ws_gc-4))]
            bins[0] = max(10, min(bins[0], self.df['cov'].max()))

        from biokit import Hist2D
        h2 = Hist2D(data)

        try:
            h2.plot(bins=bins, xlabel="Per-base coverage",
                ylabel=r'GC content (%)' ,
                Nlevels=Nlevels, contour=contour, norm=norm,
                fontsize=fontsize, **kargs)
        except:
            h2.plot(bins=bins, xlabel="Per-base coverage",
                ylabel=r'GC content (%)' ,
                Nlevels=Nlevels, contour=False, norm=norm,
                fontsize=fontsize, **kargs)
        pylab.ylim([ymin,ymax])
        corr = self.get_gc_correlation()
        return corr

    def get_gc_correlation(self):
        return self.df[['cov', 'gc']].corr().iloc[0,1]

    def get_max_gc_correlation(self, reference):
        pylab.clf()
        def func(params):
            ws = int(round(params[0]))
            if ws < 10:
                return 0
            self.bed.compute_gc_content(reference, ws)
            corr = self.get_gc_correlation()
            print(ws, corr)
            pylab.plot(ws, corr, 'o')
            return corr

        print("X, correlation\n")
        from scipy.optimize import fmin
        res = fmin(func, 200, xtol=1) # guess is 200
        pylab.xlabel("GC window size")
        pylab.ylabel("Correlation")
        return res

    def get_stats(self):
        data = self.df

        stats ={
            'DOC': self.df['cov'].mean(),
            'std': self.df['cov'].std(),
            'median': self.df['cov'].median(),
            'BOC': sum(self.df['cov'] > 0) / float(len(self.df)) }
        stats['CV'] = stats['std'] /  stats['DOC']
        stats['MAD'] = np.median(abs(data['cov'].median() - data['cov']).dropna())

        if 'scale' in self.df.columns:
            MAD = np.median(abs(data['scale'].median() - data['scale']).dropna())
            stats['MAD_normed'] = MAD
        if 'gc' in self.df.columns:
            stats['GC'] = self.df['gc'].mean() * 100

        return stats


class FilteredGenomeCov(object):
    """Class used within :class:`ChromosomeCov` to select a subset of the
    original GenomeCov

    :target: developers only
    """
    _feature_wanted = {"CDS"}
    def __init__(self, df, threshold=3, feature_list=None):
        """ .. rubric:: constructor

        :param df: dataframe with filtered position used within
            :class:`GenomeCov`. Must contain the following columns:
            ["pos", "cov", "rm", "zscore"]
        :param int threshold: size 

        """
        region_list = self._merge_region(df, threshold=threshold)
        if feature_list:
            region_list = self._add_annotation(region_list, feature_list)
        self.df = self._dict_to_df(region_list, feature_list)

    def __str__(self):
        return self.df.__str__()

    def __len__(self):
        return self.df.__len__()

    def _merge_row(self, df, start, stop):
        chrom = df["chr"][start]
        cov = np.mean(df["cov"].loc[start:stop])
        max_cov = np.max(df["cov"].loc[start:stop])
        rm = np.mean(df["rm"].loc[start:stop])
        zscore = np.mean(df["zscore"].loc[start:stop])
        max_zscore = df["zscore"].loc[start:stop].max()
        size = stop - start + 1
        return {"chr": chrom, "start": start, "end": stop + 1, "size": size,
                "mean_cov": cov, "mean_rm": rm, "mean_zscore": zscore,
                "max_zscore":max_zscore, "max_cov":max_cov}

    def _merge_region(self, df, threshold, zscore_label="zscore"):
        """Merge position side by side of a data frame.

        Uses a double threshold method.

        :param threshold: the high threshold (standard one), not the low one.

        .. todo:: to be documented
        """
        region_start = None
        region_stop = None
        start = 1
        stop = 1
        prev = 1
        # handle case where for example position n-1 have a zscore of -5 and n 
        # have a zscore of 5. It is two different regions.
        region_zscore = 0

        merge_df = []
        for pos, zscore in zip(df["pos"], df[zscore_label]):
            stop = pos
            if stop - 1 == prev and zscore * region_zscore >= 0:
                prev = stop
            else:
                if region_start:
                    merge_df.append(self._merge_row(df, region_start,
                        region_stop))
                    region_start = None
                start = stop
                prev = stop
                region_zscore = zscore
            if abs(zscore) > threshold:
                if not region_start:
                    region_start = pos
                    region_stop = pos
                else:
                    region_stop = pos

        if start < stop and region_start:
            merge_df.append(self._merge_row(df, region_start,region_stop))
        return merge_df

    def _add_annotation(self, region_list, feature_list):
        """ Add annotation from a dictionarie generated by parsers in
        sequana.tools.
        """
        region_ann = []
        # an iterator of features
        iter_feature = iter(feature_list)
        feature = next(iter_feature)
        # pass "source" feature
        while feature["type"] not in FilteredGenomeCov._feature_wanted:
            try:
                feature = next(iter_feature)
            except StopIteration:
                msg = ("Features types ({0}) are not present in the annotation "
                       "file. Please change what types you want")
                sys.exit(0)
        # merge regions and annotations
        for region in region_list:
            while feature["gene_end"] <= region["start"]:
                try:
                    feature = next(iter_feature)
                except:
                    break
            while feature["gene_start"] < region["end"]:
                # put locus_tag in gene field if gene doesn't exist
                try: 
                    feature["gene"]
                except KeyError:
                    try:
                        feature["gene"] = feature["locus_tag"]
                    except:
                        feature["gene"] = "None"
                # put note field in product if product doesn't exist
                try:
                    feature["product"]
                except KeyError:
                    try:
                        feature["product"] = feature["note"]
                    except:
                        feature["product"] = "None"
                region_ann.append(dict(region, **feature))
                try:
                    feature = next(iter_feature)
                except StopIteration:
                    break
        return region_ann

    def _dict_to_df(self, region_list, annotation):
        """ Convert dictionary as dataframe.
        """
        if annotation:
            colnames = ["chr", "start", "end", "size", "mean_cov", "mean_rm", 
                    "mean_zscore", "gene_start", "gene_end", "type", "gene", 
                    "strand", "product"]
        else:
            colnames = ["chr", "start", "end", "size", "mean_cov", "mean_rm",
                    "mean_zscore"]
        merge_df = pd.DataFrame(region_list, columns=colnames)
        int_column = ["start", "end", "size"]
        merge_df[int_column] = merge_df[int_column].astype(int)
        if annotation:
            # maybe let the user set what he wants
            return merge_df.loc[merge_df["type"].isin(
                FilteredGenomeCov._feature_wanted)]
        return merge_df

    def get_low_roi(self):
        return self.df.loc[self.df["mean_zscore"] < 0]

    def get_high_roi(self):
        return self.df.loc[self.df["mean_zscore"] > 0]

