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
import pandas as pd
from bioservices import EUtils, easyXML
import collections
from sequana.databases import ENADownload
from easydev import DevTools
import os
import ftplib
import subprocess
import sys
import glob
from easydev import execute, TempFile
import pylab
"""


kraken-build --shrink 100000 --db viruses_only --new-db test2 --shrink-block-offset 1 --kmer-len 15 --minimizer-len 8

The taxonomy files is still very large. One could filter it.


"""

class Krona(collections.Counter):
    def __init__(self, filename):
        super(Krona, self).__init__()
        self.filename = filename
        self._read()

    def _read(self, filename=None):
        if filename:
            self.filename = filename
        with open(self.filename, "r") as fin:
            for line in fin.readlines():
                count, name = line.split("\t", 1)
                self[name] += count

    def __add__(self, krona):
        self+=krona



class KrakenContaminant(object):
    """


    Run a kraken analysis. You will end up with a file e.g. kraken.out

    You could use kraken-translate but then you need extra parsing to convert
    into a Krona-compatible file. Here, we take the output from kraken and
    directly transform to the krona-compatible file.

    ::

        k = KrakenContamimant()
        k.kraken_to_krona()

    Then in a shell, ::



    """

    def __init__(self, filename="kraken.out", verbose=True):
        self.filename = filename
        self.verbose = verbose
        from biokit import Taxonomy
        self.tax = Taxonomy(verbose=self.verbose)

    def _get_taxonomy_eutils(self, ids):
        """

        "param list ids: list of taxon identifiers from whih c we want to get
            the lineage

        This won't be used anymore since it requires internet connection, not
        available on the IP cluster
        """
        if self.verbose:
            print('Connecting to NCBI/EUtils')
        e = EUtils()
        if self.verbose:
            print("Calling EFetch...")
        res = e.EFetch("taxonomy", ids)
        self.xml = easyXML(res)

        # There is one lineage per taxon so we can use the findAll
        lineage = [x.text for x in self.xml.findAll("lineage")]


        # We now want to get the scientific name for each taxon. Note, however
        # that there are several taxon children-tags within a taxon each having
        # a scientific name, so we first use getchildren() to make sure to loop
        # over a the children only
        scnames = []
        for taxon in self.xml.root.getchildren():
            name = taxon.findall("ScientificName")[0].text
            scnames.append(name)

        return lineage, scnames

    def get_taxonomy_biokit(self, ids):
        if self.verbose:
            print('Retrieving taxon using biokit.Taxonomy')

        # filter the lineage to keep only information from one of the main rank
        # that is superkingdom, kingdom, phylum, class, order, family, genus and
        # species
        ranks = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
        lineage = [self.tax.get_lineage_and_rank(x) for x in ids]

        # Now, we filter each lineage to keep only relevant ranks
        # We drop the 'no rank' and create a dictionary
        # Not nice but works for now
        results = []
        for i, this in enumerate(lineage):
            default = dict.fromkeys(ranks, ' ')
            for entry in this:
                if entry[1] in ranks:
                    default[entry[1]] = entry[0]
                elif entry[1] == "superkingdom":
                    default["kingdom"] = entry[0]
            # Scientific name is the last entry tagged as no_rank  following
            # species TODO (check this assumption)
            # e.g. 351680 and 151529 have same 7 ranks so to differenatiate
            # them, the scientific name should be used.
            # By default, we will take the last one. If species or genus, we
            # repeat the term
            try:
                default['name'] = this[-1][0]
            except:
                default['name'] = "unspecified"
            results.append(default)

        df = pd.DataFrame.from_records(results)
        df.index = ids
        df = df[list(ranks) + ['name']]

        return df

    def _parse_data(self):
        taxonomy = {}

        if self.verbose:
            print("Reading kraken data")
        # we select only col 0,2,3 to save memoty, which is required on very
        # large files
        self._df = pd.read_csv(self.filename, sep="\t", header=None, usecols=[0,2,3])

        # This gives the list of taxons as index and their amount
        # above, we select only columns 0,2,3  the column are still labelled
        # 0,2,3 in the df
        self._taxons = self._df.groupby(2).size()
        self._taxons.drop(0, inplace=True)
        self._taxons.sort_values(ascending=False, inplace=True)

        category = self.df.groupby(0).size()
        if 'C' in category.index:
            self.classified = category['C']
        else:
            pass
        if 'U' in category.index:
            self.unclassified = category['U']
        else:
            pass

    def _get_taxons(self):
        try:
            return self._taxons
        except:
            self._parse_data()
            return self._taxons
    taxons = property(_get_taxons)

    def _get_df(self):
        try:
            return self._df
        except:
            self._parse_data()
            return self._df
    df = property(_get_df)

    def kraken_to_krona(self, output_filename=None, mode=None, nofile=False):
        if output_filename is None:
            output_filename = self.filename + ".summary"
        taxon_to_find = list(self.taxons.index)

        # classified reads as root  (1)
        try:
            if self.verbose:
                print("Removing taxon 1 (%s values) " % self.taxons.ix[1])
                print("Found %s taxons " % len(taxon_to_find))
            taxon_to_find.pop(taxon_to_find.index(1))
        except:
            pass

        if mode != "adapters":
            df = self.get_taxonomy_biokit(taxon_to_find)
            self.lineage = [";".join(this) for this in df[df.columns[0:-1]].values]
            self.scnames = list(df['name'].values)  # do we need a cast ?
        else:
            # Let us get the known adapters and their identifiers
            from sequana.adapters import AdapterDB
            adapters = AdapterDB()
            adapters.load_all()

            self.scnames = []

            for taxon in self.taxons.index:
                if str(taxon) in [1, "1"]:
                    self.scnames.append('unknown')
                    continue

                if str(taxon) not in list(adapters.df.identifier):
                    self.scnames.append('unknown')
                    continue

                self.scnames.append(adapters.get_name(taxon))
            self.lineage = ["Adapters;%s"% x for x in self.scnames]

            assert len(self.lineage) == len(self.taxons)
            assert len(self.scnames) == len(self.taxons)


        with open(output_filename, "w") as fout:
            for i, this in enumerate(self.lineage):
                index = self.taxons.index[i]
                line = str(self.taxons.ix[index])+"\t"+"\t".join(this.split(';'))
                line += " " +self.scnames[i]
                fout.write(line+'\n')
            fout.write("%s\t%s" % (self.unclassified, "Unclassified"))


    def plot(self, kind="pie", cmap="copper", threshold=1, radius=0.7, **kargs):
        if kind not in ['barh', 'pie']:
            print('kind parameter: Only barh and pie are supported')
            return
        # This may have already been called but maybe not. This is not time
        # consuming, so we call it again here
        df = self.get_taxonomy_biokit(list(self.taxons.index))
        data = self.taxons.copy()
        data = data/data.sum()*100
        others = data[data<1].sum()
        data = data[data>1]

        names = df.ix[data.index]['name']
        data.index = names.values
        data.ix['others'] = others
        try:
            data.sort_values(inplace=True)
        except:
            data.sort(inplace=True)

        if kind == "pie":
            ax = data.plot(kind=kind, cmap=cmap, autopct='%1.1f%%', 
                radius=0.7, **kargs)
            pylab.ylabel(" ")
            for text in ax.texts:
                text.set_size("x-small")
        elif kind == "barh":
            ax = data.plot(kind=kind,  **kargs)
            pylab.xlabel(" percentage ")

        return data


class KrakenBuilder():
    """



    # For the viruses + subset of bacteria (salmonelle, ecoli and listeria) +
    # plasmids, this command::

        kraken-build  --rebuild -db test2 --minimizer-len 10 --max-db-size 4 --threads 4
        --kmer-len 26 --jellyfish-hash-size 500000000

    takes about 12 minutes and generates a DB of 4.2 Gb

    Downloading the taxonomy files takes about 10 minutes. This step can be
    skipped and you may just copy/paste previous downloads

    The downloads of FASTA files 


    The construction of the DB (viruses only) takes about 30-60 seconds ::

        from sequana.kraken import KrakenBuilder
        k = KrakenBuilder(dbname="viruses")
        k.download_taxonomy = True # This takes time depending on your connection
        k.run(dbs=['virus'])
 
    The directly kraken_builder should now contain a database called viruses,
    which can be used by kraken


    Kraken-build uses jellyfish. The **hash_size** parameter is the jellyfish
    hash_size parameter. If you set it to 6400M, the memort reqquired it about
    6.9bytes times 6400M that is 40Gb of memory. The default value used here
    means 3.5Gb are required. 

    The size to store the DB itself should be

    :math:

        sD + 8 (4^M)

    where **s** is about 12 bytes (used to store a kmer/taxon pair, D is the
    number of kmer in the final database, which cannot be estimated before 
    hand, and M the length minimiser parameter.

    """
    def __init__(self, dbname, download_taxon=False):
        # See databases.py module
        self.dbname = dbname
        self.enadb = ENADownload()
        self.valid_dbs = self.enadb._metadata.keys()
        self.download_taxon = download_taxon

        # mini_kraken uses minimiser-length = 13, max_db =4, others=default so
        # kmer-len=31 hashsize=default
        self.params = {
            "dbname": self.dbname,
            "minimizer_len": 10,
            "max_db_size": 4,
            "threads": 4,
            "kmer_length" : 26,
            "hash_size" : 500000000
        }

        self.init()

    def init(self):
        # mkdir library
        self._devtools = DevTools()
        self._devtools.mkdir(self.dbname)
        self._devtools.mkdir(self.dbname + os.sep + "library")
        self.added_dir = self.dbname + os.sep + "library" + os.sep + "added"
        self._devtools.mkdir(self.added_dir)
        self.taxon_dir = self.dbname + os.sep + "taxonomy" 
        self._devtools.mkdir(self.taxon_dir)

    def run(self, dbs=[]):
        """Create the Custom Kraken DB

        #. download taxonomy files
        #. Load the DBs (e.g. virus)
        #. Build DB with kraken-build

        """
        if self.download_taxon: 
            self.download_taxonomy()
        else:
            # search for taxon file. If not found, error

            required = self.taxon_dir + os.sep + "gi_taxid_nucl.dmp" 

            if required not  in glob.glob(self.taxon_dir + os.sep + "*"):
                raise IOError("Taxon file not found. Please set download_taxon to True")

        self.download_dbs(dbs)
        self.build_kraken()

    def download_taxonomy(self):
        # download taxonomy
        print('Downloading taxonomy files. Takes a while depending on your connection')
        if self.download_taxonomy:
            # We could use kraken-build --download-taxonomy + a subprocess but
            # even simpler to get the file via ftp
            FTP = "ftp.ncbi.nih.gov" 
            execute("wget %s/pub/taxonomy/gi_taxid_nucl.dmp.gz --directory-prefix %s"
                % (FTP, self.taxon_dir))

            execute("wget %s/pub/taxonomy/taxdump.tar.gz --directory-prefix %s"
                % (FTP, self.taxon_dir))

            # Unzip the files
            try:
                execute('unpigz %s/gi_taxid_nucl.dmp.gz' % self.taxon_dir)
            except:
                pass # alreaqdy done
            try:
                execute('tar xvfz %s/taxdump.tar.gz -C %s' % 
                    (self.taxon_dir, self.taxon_dir))
            except:
                pass # already done

    def download_dbs(self, dbs=[]):
        print("Downloading all Fasta files for %s" % dbs)
        # Download the DBs in it
        from .databases import ENADownload
        for db in dbs:
            assert db in self.valid_dbs, "use dbs in %s" % self.valid_dbs
            dbname, output_dir = self.enadb._metadata[db]
            self.enadb.download_fasta(dbname, self.added_dir + os.sep + output_dir)

    def build_kraken(self):
        print('Building the kraken db ')

        cmd = """kraken-build  --rebuild -db %(dbname)s \
            --minimizer-len %(minimizer_len)s\
            --max-db-size %(max_db_size)s \
            --threads %(threads)s\
            --kmer-len %(kmer_length)s \
            --jellyfish-hash-size %(hash_size)s""" % self.params

        # again, kraken-build prints on stderr so we cannot use easydev.shellcmd
        execute(cmd)


class KrakenTaxon(object):
    def __init__(self, fastq, database, threads=4, output="krona.html"):

        self.database = database
        self.threads = threads
        self.output = output

        # Fastq input
        if isinstance(fastq, str):
            self.paired = False
            self.fastq = [fastq]
        elif isinstance(fastq, list):
            if len(fastq) == 2: 
                self.paired = True
            else:
                self.paired = False
            self.fastq = fastq
        else:
            raise ValueError("Expected a fastq filename or list of 2 fastq filenames")

    def run(self):

        # We need two temp file
        kraken_summary = TempFile()
        kraken_output = TempFile()

        params = {
            "database": self.database,
            "thread": self.threads,
            "file1": self.fastq[0],
            "kraken_output": kraken_output.name,
            "kraken_summary": kraken_summary.name
            }

        if self.paired:
            params["file2"] = self.fastq[1]

        command = "kraken -db %(database)s %(file1)s " 
        if self.paired:
            command += " %(file2)s --paired"
        command += " --threads %(thread)s --out %(kraken_output)s" 


        command = command % params
        from snakemake import shell
        shell(command)

        # Translate kraken output to a format understood by Krona
        krona_summary = TempFile()
        k = KrakenContaminant(kraken_output.name)
        k.kraken_to_krona(output_filename=kraken_summary.name)

        # Transform to Krona
        shell("ktImportText %s -o %s" % (kraken_summary.name, self.output))

