from snakemake import shell as  shellcmd
import os
import shutil
import glob
import sys
from optparse import OptionParser
import argparse





class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_mapping"):
        usage = """Welcome to SEQUANA - mapping standalone

            sequana_mapping --file1 R1.fastq --file2 R2.fastq --reference reference

AUTHORS: Thomas Cokelaer
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION:
        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        # options to fill the config file
        self.add_argument("--file1", dest="file1", type=str,
            default=None, required=True,
            help="""R1 fastq file (zipped) """)
        #self.add_argument("--file2", dest="file2", type=str,
        #    default=None,
        #    help="""R2 fastq file (zipped) """)
        self.add_argument("--reference", dest="reference", type=str,
            help="""reference """)


def main(args=None):

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
       options = user_options.parse_args(args[1:])


    reference = options.reference
    if options.file1 and options.file2:
        fastq = "%s %s" % (options.file1, options.file2)
    elif options.file1 and not options.file2:
        fastq = "%s" % (options.file1)
    elif options.file1 is None:
        raise ValueError("--file1 must be used")


    fastq2 = options.file2

    params = {"reference": reference, "fastq": fastq}

    print(params)
    return
    # indexing
    shellcmd("bwa index %(reference)s.fa" % params)
    cmd = "samtools faidx %(reference)s.fa" % params

    # mapping
    cmd = r"bwa mem -t 4 -R @RG\\tID:1\\tSM:1\\tPL:illumina -T 30 %(reference)s.fa %(fastq)s  "
    cmd += "| samtools view -Sbh -> %(reference)s.bam" 
    shellcmd(cmd % params)

    # sorting BAM
    shellcmd("samtools sort -o %(reference)s.sorted.bam  %(reference)s.bam" % params)


"""reference = "JB409847"
fastq = "Cherry-1_S7_L001_R1_cutadapt_trim_1.fq"

params = {"reference": reference, "fastq": fastq}

shellcmd("bwa index %(reference)s.fa" % params)
cmd = "samtools faidx %(reference)s.fa" % params
shellcmd(r"bwa mem -t 4 -R @RG\\tID:1\\tSM:1\\tPL:illumina -T 30 %(reference)s.fa %(fastq)s  | samtools view -Sbh -> %(reference)s.bam" % params)

shellcmd("samtools sort -o %(reference)s.sorted.bam  %(reference)s.bam" % params)
"""



