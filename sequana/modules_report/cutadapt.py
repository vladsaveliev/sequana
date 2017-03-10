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
from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils import config


class CutadaptModule(SequanaBaseModule):
    """ Write HTML report of coverage analysis. This class takes either a
    :class:`GenomeCov` instances or a csv file where analysis are stored.
    """
    def __init__(self):
        """
        :param input: 
        """
        super().__init__()
        self.create_report_content()
        self.create_html("cutadapt.html")

    def create_report_content(self):
        """ Generate the sections list to fill the HTML report.
        """
        self.sections = list()
        self.add_summary_section()

    def add_summary_section(self):
        """ Coverage section.
        """
        #image = self.create_embedded_png(self.chromosome.plot_coverage,
        #                              input_arg="filename")
        self.sections.append({
            "name": "Cutadapt",
            "anchor": "cutadapt",
            "content": (
                "<p>to be done ")
        })

