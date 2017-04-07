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
import os
import io

from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils import config

from sequana.lazy import pandas as pd
from sequana.lazy import pylab
from sequana import logger

from sequana.utils.datatables_js import DataTable
from sequana.pacbio import BAMPacbio


class PacbioInputBAMModule(SequanaBaseModule):
    """ Write HTML report of Pacbio input bam. 

    """
    def __init__(self, filename, output_filename):
        """
        :param input: 
        """
        super().__init__()
        # Expected input data is the cutadapt log file
        if os.path.exists(filename) is False:
            logger.error("This file {} does not exist".format(filename))

        self.data = BAMPacbio(filename)
        self.create_report_content()
        self.create_html(output_filename)

    def create_report_content(self):
        """ Generate the sections list to fill the HTML report.
        """
        self.sections = list()
        self.add_hist_snr()

    def add_hist_snr(self):
        def plotter(filename):
                pylab.ioff()
                self.data.hist_snr()
                pylab.savefig(filename)
        html = self.create_embedded_png(plotter, "filename",
              style='width:65%')

        self.sections.append({
            "name": "Histogram individual SNRs",
            "anchor": "histogram",
            "content":  "youpi" + html
        })

