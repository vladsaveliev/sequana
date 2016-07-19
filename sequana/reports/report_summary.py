import os

from sequana.report_main import BaseReport
from easydev import DevTools

class SequanaSummary(BaseReport):
    """
    Used by the pipelines to create a summary based on the content of the
    directory. Also used by the standalone application, in which case
    config and pipeline files are not required.


    """
    def __init__(self,  directory="report", output_filename="summary.html",
                    configfile="report/config.yaml", snakefile=None, workdir=".", **kargs):

        super(SequanaSummary, self).__init__(
            jinja_filename="summary.html",
            directory=directory,
            output_filename=output_filename,
            **kargs)

        self.env.loader.searchpath.extend([workdir])

        self.workdir = workdir
        self.devtools = DevTools()

        self.title = "Summary Report"
        self.jinja['title'] = "Summary report"
        #self.jinja['description'] = "<b>Description:</b>"

        # ============================================ Add the config and pipeline files
        from sequana.snaketools import SequanaConfig
        self.config = SequanaConfig(configfile)
        self.jinja['project'] = self.config.PROJECT
        if self.config.paired:
            self.jinja['type'] = "Paired-end"
        else:
            self.jinja['type'] = "Single-end"

        # Links to the datasets
        html = ""
        for link, filename in zip(self.config.DATASET, self.config.BASENAME):
            html += '<li><a href="file:///%s">%s</a></li>\n' % (link, filename)
            html += "</ul>"
            self.jinja['dataset'] = html

        # The base has a navigation, that we do not want
        self.jinja['nav_off'] = 'True'

        if snakefile: self.read_snakefile(snakefile)

        # This is a string representation of the config file
        if configfile: self.read_configfile(configfile)


        # include whatever is relevant
        self.include_kraken()
        self.include_phix()
        self.include_sample_stats()
        self.include_adapters_stats()
        self.include_details()

        # this is a dictionary usable within JINJA templates, which may have
        # been chnaged by methods above
        self.jinja['cfg'] = self.config.config

    def include_details(self):
        self.jinja['snakemake_stats'] = "snakemake_stats.png"

    def include_adapters_stats(self):
        filename="fastq_stats__cutadapt/temp_fastq_stats__cutadapt.html"
        try:
            self.jinja['cutadapt_stats'] = open(filename, "r").read()
        except:
            self.config.config['adapter_removal']['do'] = False

    def include_sample_stats(self):
        filename="fastq_stats__samples/temp_fastq_stats__samples.html"
        self.jinja['sample_stats'] = open(filename, "r").read()

        # find an image
        import glob
        filenames = glob.glob("report/images/%s_*_R1_*boxplot*png" % self.config.PROJECT)
        if len(filenames):
            self.jinja['sample_image'] = filenames[0]

    def include_kraken(self, filename='kraken/kraken.png'):
        if os.path.exists(filename):
            self.devtools.shellcmd("cp %s report/images/kraken_pie.png" % filename)
            self.jinja['kraken_pie'] = "images/kraken_pie.png"
        else:
            self.config.config['kraken']['do'] = False

    def include_phix(self):
        filename="bwa_bam_to_fastq/bwa_mem_stats.json"
        if os.path.exists(filename):
            from sequana import tools
            stats = tools.StatsBAM2Mapped(filename)
            self.jinja['phix_section'] = stats.to_html()

        # include html pages with some stats
        filename = "fastq_stats__phix/temp_fastq_stats__phix.html"
        try:
            self.jinja['phix_stats'] = open(filename, "r").read()
        except:
            pass

    def parse(self):
        pass

