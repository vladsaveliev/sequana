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
""" Tools to launch snpEff."""
import re
import sys
import os
import shutil
import subprocess as sp
from collections import OrderedDict

from sequana.resources import snpeff
from sequana import FastA


class SnpEff(object):
    """ Python wrapper to set and launch snpEff from a genbank file. It is not
    easy to use a custom genbank file with snpEff.

    Example:
    
    ::

        snpeff = SnpEff('file.gbk')
        snpeff.launch_snpeff('variants.vcf', 'variant.ann.vcf')
    """
    extension = {"genbank": ".gbk", "gff": ".gff", "gtf": ".gtf"}
    def __init__(self, reference, file_format="", log=None):
        """.. rubric:: Constructor

        :param reference: annotation reference.
        :param file_format: format of your file. ('only genbank actually')
        :param log: log file
        """
        self.log_file = log
        if log is not None:
            if os.path.isfile(log):
                os.remove(log)
        self.reference = reference
        self.ref_name = os.path.basename(reference).split('.')[0]
        # Check if snpEff.config is present
        if not os.path.exists("snpEff.config"):
            self._get_snpeff_config()
        if not file_format:
            try:
                self._check_format()
            except FileNotFoundError:
                pass
        else:
            self.file_format = file_format
        # Check if reference is a file
        if os.path.exists(reference):
            if not os.path.exists("data" + os.sep + self.ref_name + os.sep +
                        "snpEffectPredictor.bin"):
                # Build snpEff predictor
                self._add_custom_db()
            else:
                # Add db in config if it was removed
                self._add_db_in_config()
        # Check if reference is present in snpEff database
        elif self._check_database(self.ref_name):
            if not os.path.exists("data" + os.sep + self.ref_name):
                # Download file
                snpeff_dl = sp.Popen(["snpEff", "download", self.ref_name])
                snpeff_dl.wait()
        # If reference is nowhere
        else:
            print("The file " + self.ref_name + " is not present in the "
                  "directory.\n")
            print("And your reference is not present in the database of "
                  "sequana. If you are sure that the reference is present in "
                  "the last update of the snpEff database. Please, import your "
                  "snpEff.config.\n")

    def _check_format(self):
        """ Check the format of your file.
        """
        # set regex for gff and gtf files 
        with open(self.reference, "r") as fp:
            first_line = fp.readline()
            if first_line.startswith('LOCUS'):
                self.file_format = "genbank"
            elif re.search('gff-version', first_line):
                self.file_format = "gff"
            elif first_line.startswith('#!'):
                self.file_format = "gtf"
            else:
                print("The format can not be determined, please relaunch " 
                      "the script with the file_format argument")
                sys.exit(1)

    def _check_database(self, reference):
        """ Check if your genbank is already added.
        """
        proc_db = sp.Popen(["snpEff", "databases"], stdout=sp.PIPE)
        snpeff_db = {line.split()[0] for line in proc_db.stdout}
        if reference.encode("utf-8") in snpeff_db:
            return True
        return False
    
    def _get_snpeff_config(self):
        """ Copy and unzip the snpEff.config file.
        """
        from sequana import sequana_data
        CONFIG = sequana_data("snpEff.config.gz", "snpeff")
        shutil.copyfile(CONFIG, "./snpEff.config.gz")
        gunzip_proc = sp.Popen(["gunzip", "snpEff.config.gz"])
        gunzip_proc.wait()

    def _add_custom_db(self):
        """ Add your custom file in the local snpEff database.
        """
        # create directory and copy annotation file
        genome_dir = "data" + os.sep + self.ref_name + os.sep
        try:
            os.makedirs(genome_dir)
        except FileExistsError:
            pass
        shutil.copyfile(self.reference, genome_dir + "genes" + 
                        SnpEff.extension[self.file_format])

        # add new annotation file in config file
        self._add_db_in_config()
       
        snpeff_build_line = ["snpEff", "build", "-" + self.file_format,
                             '-v', self.ref_name]
        if self.log_file:
            with open(self.log_file, "ab") as fl:
                snp_build = sp.Popen(snpeff_build_line, stderr=fl, stdout=fl)
        else:
            snp_build = sp.Popen(snpeff_build_line)
        snp_build.wait()
        rc = snp_build.returncode
        if rc != 0:
            print("snpEff build return a non-zero code")
            sys.exit(rc)

    def _add_db_in_config(self):
        """ Add new annotation at the end of snpEff.config file.
        """
        if not self._check_database(self.ref_name):
            with open("snpEff.config", "a") as fp:
                fp.write(self.ref_name + ".genome : " + self.ref_name)

    def launch_snpeff(self, vcf_filename, output, html_output=None,
                      options=""):
        """ Launch snpEff with the custom genbank file.
        
        :param str vcf_filename: input VCF filename.
        :param str output: output VCF filename.
        :param str html_output: filename of the HTML creates by snpEff.
        :param str options: any options recognised by snpEff.
        """
        # Create command line for Popen
        args_ann = ["snpEff", "-formatEff"]
        if html_output is not None:
            args_ann += ["-s", html_output]
        args_ann += [options, self.ref_name, '-v', vcf_filename]

        # Launch snpEff
        if self.log_file:
            with open(self.log_file, "ab") as fl, open(output, "wb") as fp:
                proc_ann = sp.Popen(args_ann, stdout=fp, stderr=fl)
                proc_ann.wait()
        else:
            with open(output, "wb") as fp:
                proc_ann = sp.Popen(args_ann, stdout=fp)
                proc_ann.wait()

    def _get_seq_ids(self):
        # genbank case
        if self.file_format == "genbank":
            regex = re.compile('^LOCUS\s+([\w\.\-]+)')
            chrom_regex = re.compile('\\chromosome="([\w\.\-]+)"')
            with open(self.reference, "r") as fp:
                line = fp.readline()
                seq = regex.findall(line)
                for line in fp:
                    if line.strip().startswith(('gene', 'CDS',)):
                        break
                    chrom = chrom_regex.search(line)
                    if chrom:
                         seq = [chrom.group(1)]
                         regex = chrom_regex
                seq += [regex.search(line).group(1) for line in fp 
                       if regex.search(line)]
            return seq
        # gff/gtf case
        else:
            regex = re.compile("^([^\s#]+)[ \t\v]")
            with open(self.reference, "r") as fp:
                seq = [regex.seach(line).group(1) for line in fp 
                        if regex.search(line)]
            return list(OrderedDict.fromkeys(seq))

    def add_locus_in_fasta(self, fasta, output_file):
        """ Add locus of annotation file in description line of fasta file.

        :param str fasta: input fasta file where you want to add locus.
        :param str output_file: output file.
        """
        fasta_record = FastA(fasta)
        ids_list = self._get_seq_ids()

        # check if both files have same number of contigs
        if len(fasta_record) != len(ids_list):
            print("fasta and annotation files don't have the same number of "
                  "contigs.")
            sys.exit(1)

        # check if directory exist
        output_dir = os.path.dirname(output_file)
        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        except FileNotFoundError:
            pass

        if fasta_record.names[0] == ids_list[0]:
            print("Files have same sequence id.")
            if os.path.isfile(output_file): 
                os.remove(output_file)
            os.symlink(os.path.realpath(fasta), output_file)
            return

        with open(output_file, "w") as fp:
            # write fasta with seqid of annotation file
            for n in range(len(fasta_record)):
                seq_id = ">{0} {1}\n".format(ids_list[n], fasta_record.names[n])
                seq = fasta_record.sequences[n]
                sequence = "\n".join([seq[i:min(i+80, len(seq))]
                    for i in range(0, len(seq), 80)]) + "\n"
                contigs = seq_id + sequence
                fp.write(contigs)


def download_fasta_and_genbank(identifier, tag, genbank=True, fasta=True):
    """

    :param identifier: valid identifier to retrieve from NCBI (genbank) and 
        ENA (fasta)
    :param tag: name of the filename for the genbank and fasta files.
    """
    if genbank:
        from bioservices import EUtils
        eu = EUtils()
        data = eu.EFetch(db="nuccore",id=identifier, rettype="gbwithparts",
            retmode="text")
        if isinstance(data, int) and data == 400:
            raise ValueError("%s not found on NCBI")
        else:
            with open("%s.gbk" %  tag, "w") as fout:
                fout.write(data.decode())

    if fasta:
        from bioservices import ENA
        ena = ENA()
        data = ena.get_data(identifier, 'fasta')
        if isinstance(data, int) and data == 400:
            raise ValueError("%s not found on NCBI")
        else:
            with open("%s.fa" % tag, "w") as fout:
                fout.write(data.decode())
