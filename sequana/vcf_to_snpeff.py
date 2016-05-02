""" Tools to launch snpEff.

"""

# Import -----------------------------------------------------------------------

import subprocess as sp
import re
import sys
import os
import shutil
from sequana.resources import snpeff
from sequana import sequana_data

# Global -----------------------------------------------------------------------

CONFIG = sequana_data("snpeff/snpEff.config.gz")

# Class ------------------------------------------------------------------------

class Vcf_to_snpeff(object):
    """ Python wrapper to launch snpEff.

    """
    def __init__(self, vcf_filename, reference, file_format="auto"):
        """

        :param vcf_filename: the input vcf file.
        :param reference: annotation reference.
        :param file_format: format of your file. ('-genbank'/'-gff3'/'-gtf22')
        """
        self.vcf_filename = vcf_filename
        self.reference = reference
        self.file_format = file_format
        self.ext = {"-genbank": ".gbk", "-gff3": ".gff", "-gtf22": ".gtf"}
        
        if os.path.exists("snpEff.config"):
            self.get_snpeff_config()
        
        if file_format == "auto":
            if self._check_database(reference):
                print("Everything is alright. You can launch snpEff\n")
            else:
                try:
                    self._check_format(reference)
                except IOError:
                    print("The file " + reference + " is not present in the "
                          "directory.\n")
                    print("And your reference is not present in the database "
                          "of sequana. If you are sure that the reference is "
                          "present in the last update of the snpEff database. " 
                          "Please, import your snpEff.config.\n")

    def _check_format(self, reference):
        with open(reference, "r") as fp:
            first_line = fp.readline()
            if first_line.startswith('LOCUS'):
                self.file_format = "-genbank"
            elif re.search('##gff-version +3', first_line):
                self.file_format = "-gff3"
            elif first_line.startswith('#!'):
                self.file_format = "-gtf22"
            else:
                print("The format can not be determined, please relaunch " 
                      "the script with the file_format argument")
                sys.exit(1)

    def _check_database(self, reference):
        proc_db = sp.Popen(["snpEff", "databases"], stdout=sp.PIPE)
        snpeff_db = {line.split()[0] for line in proc_db.stdout}
        if reference.encode("utf-8") in snpeff_db:
            return True
        return False
    
    def _get_snpeff_config(self):
        shutil.copyfile(CONFIG, "./snpEff.config.gz")
        gunzip_proc = sp.Popen(["gunzip", "snpEff.config.gz"])
        gunzip_proc.wait()
        
    def add_custom_db(self, stdout="build.err", stderr="build.out"):
        """ Add your custom file in the local snpEff database.

        """
        genome_dir = "data" + os.sep + self.reference + os.sep
        try:
            os.makedirs(genome_dir)
        except FileExistsError:
            pass

        shutil.copyfile(self.reference, genome_dir + "genes" + 
                self.ext[self.file_format])

        with open("snpEff.config", "a") as fp:
            print(self.reference + ".genome : " + self.reference, file=fp)
        
        with open(stdout, "wb") as out, open(stderr, "wb") as err:
            snp_build = sp.Popen(["snpEff", "build", self.file_format, 
                self.reference], stderr=err, stdout=out)
            snp_build.wait()

    def launch_snpEff(self, output, stderr="annot.err"):
        """ Launch snpEff
        
        """
        args_ann = ["snpEff", self.reference, self.vcf_filename]
        with open(output, "wb") as fp:
            proc_ann = sp.Popen(args_ann, stdout=fp)
            proc_ann.wait()
