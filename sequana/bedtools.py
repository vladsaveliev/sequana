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

from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from sequana.lazy import pylab


from sequana import running_median
from sequana.tools import gc_content
from sequana.errors import SequanaException
#from sequana.tools import genbank_features_parser


__all__ = ["GenomeCov", "ChromosomeCov", "DoubleThresholds"]


class DoubleThresholds(object):
    """Simple structure to handle the double threshold for negative and positive sides

    Used yb GenomeCov and related classes.

    ::

        dt = DoubleThresholds(-3,4,0.5,0.5)

    This means the low threshold is -3 while the high threshold is 4. The two
    following values must be between 0 and 1 and are used to define the value of
    the double threshold set to half the value of th the main threshold by default.

    Internally, the main thresholds are stored in the low and high attributes.
    The secondary thresholds are derived from the main thresholds and the
    two ratios. The ratios are named ldtr and hdtr for low double threshold ratio
    and high double threshold ration. The secondary thresholds are denoted low2 and
    high2 are are update automatically if low, high, ldtr or hdtr are changed.


    """
    def __init__(self, low=-3, high=3, ldtr=0.5, hdtr=0.5):

        assert ldtr>=0. and ldtr<=1.,\
            "ldrt parameter (low double threshold ratio) must be in [0,1]"
        assert hdtr>=0. and hdtr<=1.,\
            "hdrt parameter (high double threshold ratio) must be in [0,1]"
        assert low < 0, "low threshold must be negative"
        assert high > 0, "high threshold must be positive"

        self._ldtr = ldtr
        self._hdtr = hdtr
        self._high = high
        self._low = low

    def _get_ldtr(self):
        return self._ldtr
    def _set_ldtr(self, ldtr):
        self._ldtr = ldtr
        self._low2 = self._low * self._ldtr
    ldtr = property(_get_ldtr, _set_ldtr)

    def _get_hdtr(self):
        return self._hdtr
    def _set_hdtr(self, hdtr):
        self._hdtr = hdtr
        self._high2 = self._high * self._hdtr
    hdtr = property(_get_hdtr, _set_hdtr)

    def _get_low(self):
        return self._low
    def _set_low(self, value):
        assert value < 0.
        self._low = value
        self._low2 = self._low * self._ldtr
    low = property(_get_low,_set_low)

    def _get_high(self):
        return self._high
    def _set_high(self, value):
        assert value > 0.
        self._high = value
        self._high2 = self._high * self._ldtr
    high = property(_get_high,_set_high)

    def _get_low2(self):
        return self._low * self._ldtr
    low2 = property(_get_low2)

    def _get_high2(self):
        return self._high * self._hdtr
    high2 = property(_get_high2)

    def copy(self):
        thresholds = DoubleThresholds(self.low, self.high, 
            self.ldtr, self.hdtr)
        return thresholds

    def __str__(self):
        txt = "Low threshold: %s\n" % self.low
        txt += "High threshold: %s\n" % self.high
        txt += "double-low threshold: %s\n" % self.low2
        txt += "double-high threshold: %s" % self.high2
        return txt


class GenomeCov(object):
    """Create a list of dataframe to hold data from a BED file generated with
    samtools depth.


    This class can be used to plot the coverage resulting from a mapping, which
    is stored in BED format. The BED file may contain several chromosomes.
    There are handled independently and accessible as a list of
    :class:`ChromosomeCov` instances.

    Example:

    .. plot::
        :include-source:

        from sequana import GenomeCov, sequana_data

        filename = sequana_data('JB409847.bed')
        reference = sequana_data("JB409847.fasta")

        gencov = GenomeCov(filename)
        gencov.compute_gc_content(reference)

        gencov = GenomeCov(filename)
        for chrom in gencov:
            chrom.running_median(n=3001, circular=True)
            chrom.compute_zscore()
            chrom.plot_coverage()
        gencov[0].plot_coverage()

    Results are stored in a list of :class:`ChromosomeCov` named 
    :attr:`chr_list`.

    """

    def __init__(self, input_filename=None,
        low_threshold=-3, high_threshold=3, ldtr=0.5, hdtr=0.5):
        """.. rubric:: constructor

        :param str input_filename: the input data with results of a bedtools
            genomecov run. This is just a 3-column file. The first column is a
            string (chromosome), second column is the base postion and third
            is the coverage.
        :param float low_threshold: threshold used to identify under-covered
            genomic region of interest (ROI). Must be negative
        :param float high_threshold: threshold used to identify over-covered
            genomic region of interest (ROI). Must be positive
        :param float ldtr: fraction of the low_threshold to be used to define
            the intermediate threshold in the double threshold method. Must be
            between 0 and 1.
        :param float rdtr: fraction of the low_threshold to be used to define
            the intermediate threshold in the double threshold method. Must be
            between 0 and 1.

        """
        self.thresholds = DoubleThresholds(low_threshold, high_threshold,
            ldtr, hdtr)

        df = pd.read_table(input_filename, header=None)
        try:
            df = df.rename(columns={0: "chr", 1: "pos", 2: "cov", 3: "mapq0"})
            df = df.set_index("chr", drop=False)
            # Create a list of ChromosomeCov for each chromosome present in the
            # bedtools.
            self.chr_list = [ChromosomeCov(df.loc[key], self.thresholds) for key in
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
            if chrom.chrom_name in gc_dict.keys():
                chrom.df["gc"] = gc_dict[chrom.chrom_name]
                chrom._ws_gc = window_size
            else:
                msg = ("The chromosome (or contig) %s in your"
                       " BED/BAM file was not found in the reference provided."
                       " Make sure your input reference file is the same"
                       " as the one used to perform the mapping or just"
                       " remove the --reference parameter.")
                raise SequanaException(msg % chrom.chrom_name)

    def get_stats(self, output="json"):
        """Return basic statistics for each chromosome

        :return: dictionary with chromosome names as keys
            and statistics as values.

        .. seealso:: :class:`ChromosomeCov`.
        """
        stats = {}
        for chrom in self.chr_list:
            stats[chrom.chrom_name] = chrom.get_stats(output=output)
        return stats

    def hist(self, logx=True, logy=True, fignum=1, N=20, lw=2):
        for chrom in self.chr_list:
            chrom.plot_hist_coverage(logx=logx, logy=logy, fignum=fignum, N=N,
                histtype='step', hold=True, lw=lw)
            pylab.legend()


class ChromosomeCov(object):
    """Factory to manipulate coverage and extract region of interests.

    Example:

    .. plot::
        :include-source:

        from sequana import GenomeCov, sequana_data
        filename = sequana_data("virus.bed")

        gencov = GenomeCov(filename)

        chrcov = gencov[0]
        chrcov.running_median(n=3001)
        chrcov.compute_zscore()
        chrcov.plot_coverage()

        df = chrcov.get_roi().get_high_roi()

    The *df* variable contains a dataframe with high region of interests (over
    covered)


    .. seealso:: sequana_coverage standalone application
    """
    count = 0

    def __init__(self, df=None, thresholds=None):
        """.. rubric:: constructor

        :param df: dataframe with position for a chromosome used within
            :class:`GenomeCov`. Must contain the following columns:
            ["chr", "pos", "cov"]
        :param thresholds: a data structure :class:`DoubleThresholds` that holds
            the double threshold values.

        """
        # Chromosome position becomes the index
        ChromosomeCov.count += 1
        self.chrom_index = ChromosomeCov.count
        self.df = df.set_index("pos", drop=False)
        self.chrom_name = df["chr"].iloc[0]

        try:
            self.thresholds = thresholds.copy()
        except:
            self.thresholds = DoubleThresholds()


    def __str__(self):
        stats = self.get_stats(output="dataframe")
        stats.set_index("name", inplace=True)
        def _getter(data, key):
            return data.ix[key].Value

        txt = "\nGenome length: %s" % int(len(self.df))
        txt += "\nSequencing depth (DOC): %8.2f " % _getter(stats,'DOC')
        txt += "\nSequencing depth (median): %8.2f " % _getter(stats, 'Median')
        txt += "\nBreadth of coverage (BOC) (percent): %.2f " % _getter(stats,'BOC')
        txt += "\nGenome coverage standard deviation : %8.2f " % _getter(stats,'STD')
        txt += "\nGenome coverage coefficient variation : %8.2f " % _getter(stats,'CV')
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
        """Compute moving average of the genome coverage

        :param n: window's size. Must be odd
        :param bool circular: is the chromosome circular or not

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
        """Compute running median of genome coverage

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
        """Return the coefficient variation

        The coefficient of variation (CV) is defined as sigma / mu

        To get percentage, you must multiply by 100.

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

    def compute_zscore(self, k=2, step=10, use_em=True, verbose=True):
        """ Compute zscore of coverage and normalized coverage.

        :param int k: Number gaussian predicted in mixture (default = 2)
        :param int step: (default = 10). This parameter is used to speed
            up computation and is ignored if the length of the coverage/sequence
            is below 100,000

        Store the results in the :attr:`df` attribute (dataframe) with a
        column named *zscore*.

        .. note:: needs to call :meth:`running_median` before hand.

        """
        # here for lazy import 
        from biokit.stats import mixture
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

        # Naive checking that the 
        if k == 2:
            mus = self.gaussians['mus']
            sigmas = self.gaussians["sigmas"]

            index0 = mus.index(self.best_gaussian["mu"])
            if index0 == 0:
                mu1 = mus[1]
                s0 = sigmas[0]
                mu0 = mus[0]
            else:
                mu1 = mus[0]
                s0 = sigmas[1]
                mu0 = mus[1]
            if abs(mu0-mu1) < s0:
                if verbose:
                    print("Warning: k=2 but note that |mu0-mu1| < sigma0. ",
                        "k=1 could be a better choice") 

    def get_centralness(self):
        r"""Proportion of central (normal) genome coverage

        This is 1 - (number of non normal data) / (total length)

        .. note:: depends on the thresholds attribute being used.
        .. note:: depends slightly on :math:`W` the running median window
        """
        filtered = self.get_roi()
        return 1 - len(filtered) / float(len(self))

    def get_roi(self, features=None):
        """Keep positions with zscore outside of the thresholds range

        :param features:
        :return: a dataframe from :class:`FilteredGenomeCov`

        .. note:: depends on the :attr:`thresholds` low and high values.
        """

        try:
            second_high = self.thresholds.high2
            second_low = self.thresholds.low2
            query = "zscore > @second_high or zscore < @second_low"

            if features:
                if self.chrom_name not in features.keys():
                    features = None
                    print("Genbank name not found in list of chromosomes. " +\
    " Make sure the chromosome names in the genbank match those in the "+\
    " BED/BAM files.")
            if features:
                return FilteredGenomeCov(self.df.query(query), self.thresholds, 
                    features[self.chrom_name])
            else:
                return FilteredGenomeCov(self.df.query(query), self.thresholds)
        except KeyError:
            print("Column zscore is missing in data frame.\n"
                  "You must run compute_zscore before get low coverage.\n\n",
                  self.__doc__) 

    def plot_coverage(self, filename=None, fontsize=16, 
            rm_lw=1, rm_color="#0099cc", rm_label="Running median", 
            th_lw=1, th_color="r", th_ls="--", main_color="k", main_lw=1, 
            main_kwargs={}):

        """ Plot coverage as a function of base position.

        :param filename:

        In addition, the running median and coverage confidence corresponding to
        the lower and upper  zscore thresholds

        .. note:: uses the thresholds attribute.
        """
        # z = (X/rm - \mu ) / sigma

        high_zcov = (self.thresholds.high * self.best_gaussian["sigma"] +
                self.best_gaussian["mu"]) * self.df["rm"]
        low_zcov = (self.thresholds.low * self.best_gaussian["sigma"] +
                self.best_gaussian["mu"]) * self.df["rm"]

        pylab.clf()
        ax = pylab.gca()
        ax.set_axis_bgcolor('#eeeeee')
        pylab.xlim(0,self.df["pos"].iloc[-1])
        axes = []
        labels = []
        p1, = pylab.plot(self.df["cov"], color=main_color, label="Coverage",
                linewidth=main_lw, **main_kwargs)
        axes.append(p1)
        labels.append("Coverage")
        if rm_lw>0:
            p2, = pylab.plot(self.df["rm"],
                    color=rm_color,
                    linewidth=rm_lw,
                    label=rm_label)
            axes.append(p2)
            labels.append(rm_label)
        if th_lw>0:
            p3, = pylab.plot(high_zcov, linewidth=th_lw, color=th_color, ls=th_ls,
                label="Thresholds")
            p4, = pylab.plot(low_zcov, linewidth=th_lw, color=th_color, ls=th_ls,
                label="_nolegend_")
            axes.append(p3)
            labels.append("Thresholds")

        pylab.legend(axes, labels, loc="best")
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
            return 100
        if bins.any():
            return bins
        return 100

    def plot_hist_zscore(self, fontsize=16, filename=None, max_z=6,
            binwidth=0.5, **hist_kargs):
        """ Barplot of the zscore values

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
        """ Barplot of the normalized coverage with gaussian fitting

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

    def plot_hist_coverage(self, logx=True, logy=True, fontsize=16, N=20,
        fignum=1, hold=False, alpha=0.5, filename=None, **kw_hist):
        """


        """
        if hold is False:
            pylab.figure(fignum)
            pylab.clf()
        ax = pylab.gca()
        ax.set_axis_bgcolor('#eeeeee')

        data = self.df['cov'].dropna().values

        maxcov = data.max()
        if logx is True and logy is True:
            bins = pylab.logspace(0, pylab.log10(maxcov), N)
            pylab.hist(data, bins=bins, log=True, label=self.chrom_name,
                alpha=alpha, **kw_hist)
            pylab.semilogx()
            pylab.xlabel("Coverage (log scale)", fontsize=fontsize)
            pylab.ylabel("Count (log scale)", fontsize=fontsize)
        elif logx is False and logy is True:
            pylab.hist(data, bins=N, log=True, label=self.chrom_name,
                alpha=alpha, **kw_hist)
            pylab.xlabel("Coverage", fontsize=fontsize)
            pylab.ylabel("Count (log scale)", fontsize=fontsize)
        elif logx is True and logy is False:
            bins = pylab.logspace(0, pylab.log10(maxcov), N)
            pylab.hist(data, bins=N, label=self.chrom_name, alpha=alpha,
                **kw_hist)
            pylab.xlabel("Coverage (log scale)", fontsize=fontsize)
            pylab.ylabel("Count", fontsize=fontsize)
            pylab.semilogx()
        else:
            pylab.hist(data, bins=N, label=self.chrom_name, alpha=alpha,
                **kw_hist)
            pylab.xlabel("Coverage", fontsize=fontsize)
            pylab.ylabel("Count", fontsize=fontsize)
        pylab.grid(True)
        if filename:
            pylab.savefig(filename)

    def write_csv(self, filename, start=None, stop=None, header=True):
        """ Write CSV file of the dataframe.

        :param filename: csv output filename.
        :param header: boolean which determinate if the header is written.

        """
        try:
            labels=["pos", "cov", "mapq0"]
            self.df[labels][start:stop].to_csv(filename, header=header)
        except NameError:
            print("You must set the file name")
        except KeyError:
            self.df[labels[:-1]][start:stop].to_csv(filename, header=header)

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
        """Return the correlation between the coverage and GC content

        The GC content is the one computed in :meth:`GenomeCov.compute_gc_content`
        (default window size is 101)

        """
        return self.df[['cov', 'gc']].corr().iloc[0,1]

    def get_max_gc_correlation(self, reference):
        """Plot correlation between coverage and GC content by varying the GC window

         The GC content uses a moving window of size W. This parameter affects
         the correlation bewteen coverage and GC. This function find the *optimal*
         window length.

        """
        pylab.clf()
        corrs = []
        wss = []
        def func(params):
            ws = int(round(params[0]))
            if ws < 10:
                return 0
            self.bed.compute_gc_content(reference, ws)
            corr = self.get_gc_correlation()
            corrs.append(corr)
            wss.append(ws)
            return corr

        from scipy.optimize import fmin
        res = fmin(func, 100, xtol=1, disp=False) # guess is 200
        pylab.plot(wss, corrs, "o")
        pylab.xlabel("GC window size")
        pylab.ylabel("Correlation")
        pylab.grid()
        return res[0]

    def get_stats(self, output="json"):
        """Return basic stats about the coverage data"""
        data = self.df

        stats ={
            'DOC': self.df['cov'].mean(),
            'STD': self.df['cov'].std(),
            'Median': self.df['cov'].median(),
            'BOC': 100 * sum(self.df['cov'] > 0) / float(len(self.df)) }
        stats['CV'] = stats['STD'] /  stats['DOC']
        stats['MAD'] = np.median(abs(data['cov'].median() - data['cov']).dropna())

        names = ['BOC','CV', 'DOC', 'MAD', 'Median', "STD"]
        descriptions = [("breadth of coverage: the proportion (in %s) of the "
            "genome covered by at least one read. "),
            ("the coefficient of variation"),
            ("the sequencing depth (Depth of Coverage), that is the average of"
            "the genome coverage"),
            ("median of the absolute median deviation median(|X-median(X)|)"),
            ("Median of the coverage"),
            ("standard deviation")
        ]

        if 'gc' in self.df.columns:
            stats['GC'] = self.df['gc'].mean() * 100
            names.append('GC')
            descriptions.append("GC content in %")

        df = pd.DataFrame({
            "name":names, 
            "Value":[stats[x] for x in names],
            "Description": descriptions})

        if output == "json":
            return df.to_json()
        else:
            return df


class FilteredGenomeCov(object):
    """Class used within :class:`ChromosomeCov` to select a subset of the
    original GenomeCov

    :target: developers only
    """
    _feature_not_wanted = {"gene", "regulatory", "source"}
    def __init__(self, df, threshold, feature_list=None):
        """ .. rubric:: constructor

        :param df: dataframe with filtered position used within
            :class:`GenomeCov`. Must contain the following columns:
            ["pos", "cov", "rm", "zscore"]
        :param int threshold: a :class:`~sequana.bedtools.DoubleThresholds` instance.

        """
        if isinstance(feature_list, list) and len(feature_list) == 0:
            feature_list = None
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
        if zscore >= 0 :
            max_zscore = df["zscore"].loc[start:stop].max()
        else:
            max_zscore = df["zscore"].loc[start:stop].min()
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

            if zscore >0 and  zscore> threshold.high:
                if not region_start:
                    region_start = pos
                    region_stop = pos
                else:
                    region_stop = pos
            elif zscore <0 and zscore<threshold.low:
                if not region_start:
                    region_start = pos
                    region_stop = pos
                else:
                    region_stop = pos

        if start < stop and region_start:
            merge_df.append(self._merge_row(df, region_start,region_stop))
        return merge_df

    def _add_annotation(self, region_list, feature_list):
        """ Add annotation from a dictionary generated by parsers in
        sequana.tools.
        """
        region_ann = []
        # an iterator of features
        iter_feature = iter(feature_list)
        feature = next(iter_feature)
        # pass "source" feature
        while feature["type"] in FilteredGenomeCov._feature_not_wanted:
            try:
                feature = next(iter_feature)
            except StopIteration:
                msg = ("Features types ({0}) are not present in the annotation "
                       "file. Please change what types you want")
                return region_ann
        # merge regions and annotations
        for region in region_list:
            feature_exist = False
            while feature["gene_end"] <= region["start"]:
                try:
                    feature = next(iter_feature)
                except:
                    break
            while feature["gene_start"] < region["end"]:
                # A feature exist for detected ROI
                feature_exist = True
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
                if region["start"] == 237433:
                    print(dict(region, **feature))
                region_ann.append(dict(region, **feature))
                try:
                    feature = next(iter_feature)
                except StopIteration:
                    break
            if feature_exist is False:
                region_ann.append(dict(region, **{"gene_start": None,
                                                  "gene_end": None,
                                                  "type": None,
                                                  "gene": None,
                                                  "strand": None,
                                                  "product": None}))
        return region_ann

    def _dict_to_df(self, region_list, annotation):
        """ Convert dictionary as dataframe.
        """
        if annotation:
            colnames = ["chr", "start", "end", "size", "mean_cov", "mean_rm", 
                        "mean_zscore", "max_zscore", "gene_start", "gene_end",
                        "type", "gene", "strand", "product"]
        else:
            colnames = ["chr", "start", "end", "size", "mean_cov", "mean_rm",
                    "mean_zscore", "max_zscore"]
        merge_df = pd.DataFrame(region_list, columns=colnames)
        int_column = ["start", "end", "size"]
        merge_df[int_column] = merge_df[int_column].astype(int)
        if annotation:
            merge_df.rename(columns = {"gene": "gene_name"}, inplace=True)
            # maybe let the user set what he wants
            return merge_df.loc[~merge_df["type"].isin(
                FilteredGenomeCov._feature_not_wanted)]
        return merge_df

    def get_low_roi(self):
        return self.df.loc[self.df["max_zscore"] < 0]

    def get_high_roi(self):
        return self.df.loc[self.df["max_zscore"] >= 0]

