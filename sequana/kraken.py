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
from easydev import execute

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

    def __init__(self, filename="kraken.out"):
        self.filename = filename

    def get_taxonomy(self, ids):
        """

        "param list ids: list of taxon identifiers from whih c we want to get
            the lineage
        """
        print('Connecting to NCBI/EUtils')
        e = EUtils()
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

    def _parse_data(self):

        # Each entry is made of an identifier and a taxonomy
            # The taxonomy looks like:
        #     d__Viruses|o__Mononegavirales|f__Paramyxoviridae|g__Morbillivirus|s__Measles_virus
        # where we can identify the standard levels of taxonomy with rank
        # assignments from kingdom, phylum, class, order, family, genus, species that
        # can be identifier with leading character as d, p, c, o, f, g, s
        # Note that "d" stands for kingdom and superkindom is dropped
        tags = ["d", "p", "c", "o", "f", "g", "s"]

        taxonomy = {}

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

    def kraken_to_krona(self, output_filename=None, mode=None):
        if output_filename is None:
            output_filename = self.filename + ".summary"
        taxon_to_find = list(self.taxons.index)
        print("Fetching %s taxons " % len(taxon_to_find))

        if mode != "adapters":
            self.scnames = []
            self.lineage = []
            from easydev import split_into_chunks
            N = int(len(taxon_to_find)/180.) + 1

            for chunk in split_into_chunks(taxon_to_find, N):
                lineage, scnames = self.get_taxonomy(chunk)
                self.lineage.extend(lineage)
                self.scnames.extend(scnames)
        else:
            # Let us get the known adapters and their identifiers
            from sequana.adapters import AdapterDB
            adapters = AdapterDB()
            adapters.load_all()

            self.scnames = [adapters.get_name(ide) if ide not in [1, "1"]
                else "unknown" for ide in self.taxons.index]
            self.lineage = ["Adapters;%s" for x in self.scnames]

            assert len(self.lineage) == len(self.taxons)
            assert len(self.scnames) == len(self.taxons)


        with open(output_filename, "w") as fout:
            for i, this in enumerate(self.lineage):
                index = self.taxons.index[i]
                line = str(self.taxons.ix[index])+"\t"+"\t".join(this.split(';'))
                line += " " +self.scnames[i]
                fout.write(line+'\n')
            fout.write("%s\t%s" % (self.unclassified, "Unclassified"))





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

