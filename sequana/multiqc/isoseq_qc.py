#!/usr/bin/env python

""" MultiQC module to parse output from sequana"""
import os
import re

# prevent boring warning (version 1.0)
import logging
logging.captureWarnings(True)
from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table, heatmap, bargraph
logging.captureWarnings(False)

# Initialise the logger
log = logging.getLogger('multiqc.sequana/isoseq_qc')


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/isoseq_qc',    # name that appears at the top
            anchor='sequana_isoseq_qc',  # ??
            target='sequana_isoseq_qc',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="pipelines multi Summary")

        self.sequana_data = {}
        for myfile in self.find_log_files("sequana/isoseq_qc", filehandles=True):
            #print( myfile['f'] )       # File contents
            #print( myfile['s_name'] )  # Sample name (from cleaned filename)
            #print( myfile['fn'] )      # Filename
            #print( myfile['root'] )    # Directory file was in
            name = myfile['s_name']
            print(name)
            #if name.startswith("summary_"):
            #    name = name.replace("summary_", "")

            try:
                parsed_data = self.parse_logs(myfile["f"])
            except:
                print("{} could not be parsed".format(name))
                continue
            name = parsed_data["s_name"]
            self.sequana_data[name] = parsed_data

        info = "<ul>"
        for this in sorted(self.sequana_data.keys()):
            info += '<li><a href="{}/summary.html">{}</a></li>'.format(this,this,this)
        info += "</ul>"
        href="http://sequana.readthedocs.io/en/master/"
        target = "Sequana"
        mname = '<a href="{}" target="_blank">{}</a> individual report pages:'.format(href, target)
        self.intro = '<p>{} {}</p>'.format( mname, info)

        if len(self.sequana_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.sequana_data)))

        self.populate_columns()
        self.add_productivity()

    def add_productivity(self):
        from collections import OrderedDict
        data = OrderedDict()
        for name in self.sequana_data.keys():
            data[name] = {
                'P0': self.sequana_data[name]["P0"],
                'P1': self.sequana_data[name]["P1"],
                'P2': self.sequana_data[name]["P2"]
            }

        pconfig = {
            "title": "Productivity",
            "percentages": False,
            "min": 100,
            "logswitch": True,
        }
        self.add_section(
            name = 'Productivity',
            anchor = 'productivity',
            description = 'P0/P1/P2 productivities',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def parse_logs(self, log_dict):
        import json
        log_dict = json.loads(log_dict)
        data = {}
        #data['count'] = log_dict['data']["count"]
        data["P0"] = log_dict['data']['QC']["productivity"]["P0"]
        data["P1"] = log_dict['data']['QC']["productivity"]["P1"]
        data["P2"] = log_dict['data']['QC']["productivity"]["P2"]
        data["alldata"] = log_dict
        data["s_name"] = log_dict['sample_name']

        return data

    def populate_columns(self):
        from collections import OrderedDict
        headers = OrderedDict()

        for this in ["P0", "P1", "P2"]:
            if any([this in self.sequana_data[s] for s in self.sequana_data]):
               headers[this] = {
                    'title': this,
                    'description': this,
                    'min': 0,
                    'max': 100,
                    'scale': 'RdYlGn',
                    'format': '{:,.0d}',
                    'shared_key': 'count',
                }

        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)


