import pandas as pd
from bioservices import EUtils, easyXML
import collections


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
            from easydev.chunks import baskets_from
            N = int(len(taxon_to_find)/180.) + 1

            for chunk in baskets_from(taxon_to_find, N):
                lineage, scnames = self.get_taxonomy(chunk)
                self.lineage.extend(lineage)
                self.scnames.extend(scnames)
        else:
            # Let us get the known adapters and their identifiers
            from sequana.adapters import AdapterDB
            adapters = AdapterDB()
            adapters.load_all()

            self.scnames = [adapters.get_name(ide) if ide not in [1, "1"] else
"unknown" for ide in self.taxons.index]
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



    """
    def __init__(self):
        # See databases.py module 
        self.valid_dbs = ["virus", "bacteria", "plasmids"]


    def builder(self, dbs=[]):
        # Create a local direcory

        # Download the DBs in it

        # launch the command with proper options
        pass
