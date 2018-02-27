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
log = logging.getLogger('multiqc.sequana/isoseq')


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):
        print("------------------------------------")
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/isoseq',    # name that appears at the top
            anchor='sequana_isoseq',  # ??
            target='sequana_isoseq',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="pipelines multi Summary")

        print("------------------------------------")
        self.sequana_data = {}
        files = self.find_log_files("sequana/isoseq", filehandles=True)
        print(len(list(files)))
        for myfile in self.find_log_files("sequana/isoseq"):
            #print( myfile['f'] )       # File contents
            #print( myfile['s_name'] )  # Sample name (from cleaned filename)
            #print( myfile['fn'] )      # Filename
            #print( myfile['root'] )    # Directory file was in
            name = myfile['s_name']
            print("test")
            #if name.startswith("summary_"):
            #    name = name.replace("summary_", "")

            try:
                parsed_data = self.parse_logs(myfile["f"])
            except:
                print("{} could not be parsed".format(name))
                continue
            name = parsed_data["s_name"]
            self.sequana_data[name] = parsed_data
        print("------------------------------------")

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
        self.add_ccs_reads_section()
        self.add_ccs_mean_length_section()

    def add_ccs_reads_section(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                'number_ccs_reads': self.sequana_data[name]["number_ccs_reads"]
            }

        pconfig = {
            "title": "Number of CCS reads",
            "percentages": False,
            "min": 100,
            "logswitch": True,
        }
        self.add_section(
            name = 'Number of CCS reads',
            anchor = 'number_ccs_reads',
            description = 'Number of CCS reads',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def add_ccs_mean_length_section(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                'mean': self.sequana_data[name]["mean_length"],
            }
        pconfig = {
                "title": "Mean CCS read length",
            "percentages": False,
            "min": 100,
            "logswitch": True,
        }

        self.add_section(
            name = 'Mean CCS read length',
            anchor = 'mean_ccs_read_length',
            description = 'Mean CCS length of the reads',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))

    def parse_logs(self, log_dict):
        import json
        log_dict = json.loads(log_dict)
        data = {}
        #data['count'] = log_dict['data']["count"]
        data["mean_length"] = log_dict['data']['CCS']["mean_length"]
        data["number_ccs_reads"] = log_dict['data']['CCS']["number_ccs_reads"]
        data["alldata"] = log_dict
        data["s_name"] = log_dict['sample_name']

        #data["mean_gc"] = log_dict['mean_gc']
        #data["hist_gc"] = log_dict['hist_gc']
        #data["hist_read_length"] = log_dict['hist_read_length']
        return data

    def populate_columns(self):
        headers = {}
        if any(['mean_length' in self.sequana_data[s] for s in self.sequana_data]):
            headers['mean_length'] = {
                'title': 'CCS mean read length',
                'description': 'CCS mean read length',
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.0d}',
                'shared_key': 'count',
            }

        if any(['number_ccs_reads' in self.sequana_data[s] for s in self.sequana_data]):
            headers["number_ccs_reads"] = {
                'title': 'CCS reads',
                'description': 'Number of CCS reads',
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.0d}'
            }

        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)


