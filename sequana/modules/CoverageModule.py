# coding: utf-8
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
"""Module to write coverage report"""
from sequana import bedtools
from sequana.modules import GenericModule


class CoverageModule(object):
    """ Write HTML report of coverage analysis. This class takes either a
    :class:`GenomeCov` instances or a csv file where analysis are stored.
    """
    def __init__(self, input_file, template, output_directory):
        """
        :param input: 
        """
        try:
            self.bed = bedtools.GenomeCov(input_file)
        except IOError:
            msg = ("The csv file is not present. Please, check if your"
                   " file is present.")
            raise IOError(msg)
        except TypeError:
            self.bed = input_file
        # try:
        self.create_reports(template, output_directory)
        #except:
        #    msg = ("Input must be either a csv file or a :class:`GenomeCov` "
        #           "instance.")
        #    raise TypeError(msg)

    def create_reports(self, template, output_dir):
        """ Create HTML report for each chromosome present in data.
        """
        for chrom in self.bed:
            jinja = ChromosomeCoverageModule(chrom, template, output_dir)


class ChromosomeCoverageModule(GenericModule.GenericModule):
    """ Write HTML report of coverage analysis for each chromosome. It is
    created by CoverageModule.
    """
    def __init__(self,
                 chromosome,
                 template,
                 output_directory="report/"):
        """
        """
        super().__init__(template=template,
                         output_directory=output_directory)
        self.chromosome = chromosome
        self.create_report_content()
        self.create_html(template, output_directory)

    def create_report_content(self):
        """
        """
        self.sections = list()
        
        self.coverage_plot()
        self.coverage_barplot()
        # self.submappings()
        self.basic_stats()
        self.regions_of_interest()
        # self.normalized_coverage()
        # self.zscore_distribution()

    def create_html(self, template, output_directory):
        """ Create html with Jinja2.
        """
        report_output = self.j_template.render(base_content="content.html",
                                               module=self)
        with open(output_directory + "test.html", "w") as fp:
            print(report_output, file=fp)

    def coverage_plot(self):
        """
        """
        image = self.create_embed_png(self.chromosome.plot_coverage)
        self.sections.append({
            "name": "Coverage",
            "anchor": "coverage",
            "content": (
                "<p>The following figure shows the per-base coverage along the"
                " reference genome (black line). The blue line indicates the "
                "running median. From the normalised coverage, we estimate "
                "z-scores on a per-base level. The red lines indicates the "
                "z-scores at plus or minus N standard deviations, where N is "
                "chosen by the user (default:4)</p>\n{0}".format(image))
        })

    def coverage_barplot(self):
        """
        """
        image1 = self.create_embed_png(self.chromosome.plot_hist_coverage,
                                       style="width:45%")
        image2 = self.create_embed_png(self.chromosome.plot_hist_coverage,
                                       kwargs={"logx": False},
                                       style="width:45%")
        self.sections.append({
            "name": "Coverage histogram",
            "anchor": "cov_barplot",
            "content": (
                "<p>The following figure contains the histogram of the genome "
                "coverage. The X and Y axis being in log scale.</p>\n"
                "{0}\n{1}".format(image1, image2))
        })

    def basic_stats(self):
        """
        """
        html_table = self.dataframe_to_html_table(
            self.chromosome.get_stats(output="dataframe"), {"index": False})
        self.sections.append({
            "name": "Basic stats",
            "anchor": "basic_stats",
            "content": (
                "<p>The following table gives some basic statistics about the "
                "genome coverage.</p>\n{0}".format(html_table))
        })

    def regions_of_interest(self):
        """
        """
        rois = self.chromosome.get_roi()
        low_roi = rois.get_low_roi()
        high_roi = rois.get_high_roi()
        html_low_roi = self.dataframe_to_html_table(low_roi)
        html_high_roi = self.dataframe_to_html_table(high_roi)
        roi_paragraph = (
            "Regions with a z-score {0}er than {1:.2f} and at "
            "least one base with a z-score {0}er than {2:.2f} are detected as "
            "{0} coverage region. Thus, there are {3} {0} coverage regions.")
        low_paragraph = roi_paragraph.format("low",
                                             self.chromosome.thresholds.low2,
                                             self.chromosome.thresholds.low,
                                             len(low_roi))
        high_paragraph = roi_paragraph.format("high",
                                              self.chromosome.thresholds.high2,
                                              self.chromosome.thresholds.high,
                                              len(high_roi))
        self.sections.append({
            "name": "Regions Of Interest (ROI)",
            "anchor": "roi",
            "content": (
                "<p>Running median is the median computed along the genome "
                "using a sliding window. The following tables give regions of "
                "interest detected by sequana. Here is some captions:</p>\n"
                "<ul><li>mean_cov: the average of coverage</li>\n"
                "<li>mean_rm: the average of running median</li>\n"
                "<li>mean_zscore: the average of zscore</li>\n"
                "<li>max_zscore: the higher zscore contains in the region</li>"
                "</ul>\n"
                "<h3>Low coverage region</h3>\n{0}\n{1}\n"
                "<h3>High coverage region</h3>\n{2}\n{3}\n").format(
                low_paragraph, html_low_roi, high_paragraph, html_high_roi)
        })
