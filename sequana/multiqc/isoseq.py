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

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='Sequana/isoseq',    # name that appears at the top
            anchor='sequana',  # ??
            target='sequana',  # Name show that link to the following href
            href='http://github.com/sequana/sequana/',
            info="pipelines multi Summary")

        self.sequana_data = {}
        for myfile in self.find_log_files("sequana/isoseq"):
            #print( myfile['f'] )       # File contents
            #print( myfile['s_name'] )  # Sample name (from cleaned filename)
            #print( myfile['fn'] )      # Filename
            #print( myfile['root'] )    # Directory file was in
            name = myfile['s_name']
            #print(name)
            #if name.startswith("summary_"):
            #    name = name.replace("summary_", "")

            parsed_data = self.parse_logs(myfile["f"])
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
        self.add_ccs_section()

    def add_ccs_section(self):
        data = {}
        for name in self.sequana_data.keys():
            data[name] = {
                'mean': self.sequana_data[name]["mean_length"],
                'number_ccs_reads': self.sequana_data[name]["number_ccs_reads"]
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

        pconfig['title'] = "Number of CCS reads"
        self.add_section(
            name = 'Number of CCS reads',
            anchor = 'number_ccs_reads',
            description = 'Number of CCS reads',
            helptext = "",
            plot = bargraph.plot(data, None, pconfig))


    def add_hist_GC(self):
        """ Create the HTML for the FastQC GC content plot """
        data = dict()
        data_norm = dict()
        for s_name in self.sequana_data:
            try:
                X = self.sequana_data[s_name]["hist_gc"]['X']
                Y = self.sequana_data[s_name]["hist_gc"]['Y']
                Y = [y / float(sum(Y)) for y in Y]
                data[s_name] = {x:10*y for x,y in zip(X[1:], Y)}
            except KeyError:
                pass
            #else:
            #    data_norm[s_name] = dict()
            #    total = sum( [ c for c in data[s_name].values() ] )
            #    for gc, count in data[s_name].items():
            #        data_norm[s_name][gc] = (count / total) * 100

        if len(data) == 0:
            log.debug('no data for the GC content plots')
            return None

        pconfig = {
            'id': 'sequana_pacbio_per_sequence_gc_content_plot',
            'title': 'Per Sequence GC Content',
            'ylab': 'Count',
            'xlab': '% GC',
            'ymin': 0,
            #'ymax': 0.2,
            'xmax': 100,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}',
            #'colors': self.get_status_cols('per_sequence_gc_content'),
            'data_labels': [
                {'name': 'Percentages', 'ylab': 'Percentage'},
                {'name': 'Counts', 'ylab': 'PDF'}
            ]
        }

        self.add_section (
            name = 'Per Sequence GC Content',
            anchor = 'fastqc_per_sequence_gc_content',
            description = "GC content (normalised)",
            #plot = linegraph.plot([data_norm, data], pconfig))
            plot = linegraph.plot(data, pconfig))

    def add_hist_length(self):
        """ Create the HTML for the FastQC GC content plot """
        data = dict()
        data_norm = dict()
        for s_name in self.sequana_data:
            X = self.sequana_data[s_name]["hist_read_length"]['X']
            Y = self.sequana_data[s_name]["hist_read_length"]['Y']
            #Y = [y / sum(Y) for y in Y]
            data[s_name] = {x:y for x,y in zip(X[1:], Y)}
            try:
                X = self.sequana_data[s_name]["hist_read_length"]['X']
                Y = self.sequana_data[s_name]["hist_read_length"]['Y']
                #Y = [y / sum(Y) for y in Y]
                data[s_name] = {x:y for x,y in zip(X[1:], Y)}
            except KeyError:
                pass

        if len(data) == 0:
            log.debug('no data for the read length plots')
            return None

        pconfig = {
            'id': 'sequana_pacbio_hist_length',
            'title': 'Per Sequence GC Content',
            'ylab': '#',
            'xlab': 'Length',
            'ymin': 0,
            'xmax': 50000,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}Length</b>: {point.y}',
            #'colors': self.get_status_cols('per_sequence_gc_content'),
            #'data_labels': [
            #    {'name': 'Percentages', 'ylab': 'Percentage'},
            #    {'name': 'length', 'ylab': '#'}
            #]
        }

        self.add_section (
            name = 'Read length histograms',
            anchor = 'fastqc_per_sequence_gc_content',
            description = "GC content (normalised)",
            plot = linegraph.plot(data, pconfig))

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
                'max': 100,
                'min': 0,
                'scale': 'RdYlGn',
                'format': '{:,.2f}'
            }

        print(headers)
        if len(headers.keys()):
            self.general_stats_addcols(self.sequana_data, headers)


