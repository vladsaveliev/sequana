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
import os

from easydev import DevTools
from easydev import execute, TempFile

import pandas as pd

import pylab

from sequana.misc import wget
from sequana import sequana_config_path


__all__ = ['KrakenResults', "KrakenPipeline", "KrakenAnalysis"]


class KrakenResults(object):
    """Translate Kraken results into a Krona-compatible file


    If you run a kraken analysis with :class:`KrakenAnalysis`, you will end up
    with a file e.g. named kraken.out (by default).

    You could use kraken-translate but then you need extra parsing to convert
    into a Krona-compatible file. Here, we take the output from kraken and
    directly transform it to a krona-compatible file.

    ::

        k = KrakenResults("kraken.out")
        k.kraken_to_krona()

    .. note:: This takes care of fetching taxons and the corresponding lineages
        from online web services.

    """
    def __init__(self, filename="kraken.out", verbose=True):
        """

        :param filename: the input from KrakenAnalysis class
        :param verbose:

        """
        self.filename = filename
        self.verbose = verbose

        on_rtd = os.environ.get("READTHEDOCS", None) == "True"

        if on_rtd is False:
            from biokit import Taxonomy
            self.tax = Taxonomy(verbose=self.verbose)
        else:
            class Taxonomy(object):
                from sequana import sequana_data # must be local
                df = pd.read_csv(sequana_data("test_taxon_rtd.csv"), 
                    index_col=0)
                def get_lineage_and_rank(self, x):
                    # Note that we add the name as well here
                    ranks = ['kingdom', 'phylum', 'class', 'order', 
                            'family', 'genus', 'species', 'name']
                    return [(self.df.ix[x][rank], rank) for rank in ranks]
            self.tax = Taxonomy()

        # This initialise the data
        self._parse_data()

        self._data_created = False

    def get_taxonomy_biokit(self, ids):
        """Retrieve taxons given a list of taxons

        :param list ids: list of taxons as strings or integers. Could also
            be a single string or a single integer
        :return: a dataframe
        .. note:: the first call first loads all taxons in memory and takes a
            few seconds but subsequent calls are much faster
        """
        if self.verbose:
            print('Retrieving taxon using biokit.Taxonomy')

        if isinstance(ids, list) is False:
            ids = [ids]

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

        df.index = df.index.astype(int)

        return df

    def _parse_data(self):
        taxonomy = {}

        if self.verbose:
            print("Reading kraken data")
        # we select only col 0,2,3 to save memoty, which is required on very
        # large files
        self._df = pd.read_csv(self.filename, sep="\t", header=None,
                               usecols=[0,2,3])

        # This gives the list of taxons as index and their amount
        # above, we select only columns 0,2,3  the column are still labelled
        # 0,2,3 in the df
        self._taxons = self._df.groupby(2).size()
        try:
            self._taxons.drop(0, inplace=True)
        except:
            pass # 0 may not be there
        self._taxons.sort_values(ascending=False, inplace=True)

        category = self.df.groupby(0).size()
        if 'C' in category.index:
            self.classified = category['C']
        if 'U' in category.index:
            self.unclassified = category['U']

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

        if len(taxon_to_find) == 0:
            print("No reads were identified. You will need a more complete database")
            return
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

        # Now save the file
        self.output_filename = output_filename
        with open(output_filename, "w") as fout:
            for i, this in enumerate(self.lineage):
                index = self.taxons.index[i]
                line = str(self.taxons.ix[index])+"\t"+"\t".join(this.split(';'))
                line += " " +self.scnames[i]
                fout.write(line+'\n')
            try:
                fout.write("%s\t%s" % (self.unclassified, "Unclassified"))
            except:
                pass #unclassified may not exists
        self._data_created = True

    def plot(self, kind="pie", cmap="copper", threshold=1, radius=0.7, **kargs):
        """A simple non-interactive plot of taxons

        A Krona Javascript output is also available in :meth:`kraken_to_krona`

        .. plot::
            :include-source:

            from sequana import KrakenResults, sequana_data
            test_file = sequana_data("test_kraken.out", "testing")
            k = KrakenResults(test_file, verbose=False)
            df = k.plot(kind='pie')

        .. seealso:: to generate the data see :class:`KrakenPipeline`
            or the standalone application **sequana_taxonomy**.
        """
        if self._data_created == False:
            self.kraken_to_krona()

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
        
    def to_js(self, output="krona.html", onweb=False):
        if self._data_created == False:
            self.kraken_to_krona()
        execute("ktImportText %s -o %s" % (self.output_filename, output))
        if onweb is True:
            import easydev
            easydev.onweb(output)


class KrakenPipeline(object):
    """Used by the standalone application sequana_taxonomy

    This runs Kraken on a set of FastQ files, transform the results
    in a format compatible for Krona, and creates a Krona HTML report.

    ::

        from sequana import KrakenTaxon
        kt = KrakenPipeline(["R1.fastq.gz", "R2.fastq.gz"], database="krakendb")
        kt.run()
        kt.show()

    .. warning:: We do not provide Kraken database within sequana. You may
        either download a database from https://ccb.jhu.edu/software/kraken/
        or use this class to download a toy example that will
        be stored in e.g .config/sequana under Unix platforms.
        See :class:`KrakenDownload`.


    .. seealso:: We provide a standalone application of this class, which is
        called sequana_taxonomy and can be used within a command shell.

    """
    def __init__(self, fastq, database, threads=4, output="krona.html"):
        """.. rubric:: Constructor

        :param fastq: either a fastq filename or a list of 2 fastq filenames
        :param database: the path to a valid Kraken database
        :param threads: number of threads to be used by Kraken
        :param output: output filename of the Krona HTML page

        Description: internally, once Kraken has performed an analysis, reads
        are associated to a taxon (or not). We then find the correponding
        lineage and scientif names to store within a Krona formatted file.
        KtImportTex is then used to create the Krona page.

        """
        self.ka = KrakenAnalysis(fastq, database, threads)
        self.output = output

    def run(self, keep_temporary_files=False):
        """Run the analysis using Kraken and create the Krona output"""

        # Run Kraken
        self.ka.run()

        # Translate kraken output to a format understood by Krona
        kraken_summary = TempFile()
        kr = KrakenResults(self.ka.kraken_output.name)
        kr.kraken_to_krona(output_filename=kraken_summary.name)

        # Transform to Krona
        from snakemake import shell
        shell("ktImportText %s -o %s" % (kraken_summary.name, self.output))

        print(self.ka.kraken_output.name)
        print(kraken_summary.name)
        #if keep_temporary_files is False:
        #    self.ka.kraken_output.delete()
        #    kraken_summary.delete()

    def show(self):
        """Opens the filename defined in the constructor"""
        from easydev import onweb
        onweb(self.output)


class KrakenAnalysis(object):
    """Run kraken on a set of FastQ files

    In order to run a Kraken analysis, we firtst need a local database.
    We provide a Toy example. The ToyDB is downloadable as follows ( you will
    need to run the following code only once)::

        from sequana import KrakenDownload
        kd = KrakenDownload()
        kd.download_kraken_toydb()

    .. seealso:: :class:`KrakenDownload`  for more database and
        :class:`sequana.kraken_builder.KrakenBuilder` to build your own
        databases

    The path to the database is required to run the analysis. It has been
    stored in the directory ./config/sequana/kraken_toydb under Linux platforms
    The following code should be platform independent::

        import os
        from sequana import sequana_config_path
        database = sequana_config_path + os.sep + "kraken_toydb")

    Finally, we can run the analysis on the toy data set::

        from sequana import sequana_data
        data = sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
        ka = KrakenAnalysis(data, database=database)
        ka.run()

    This creates a file named *kraken.out*. It can be interpreted with
    :class:`KrakenResults`
    """
    def __init__(self, fastq, database, threads=4 ):
        """.. rubric:: Constructor

        :param fastq: either a fastq filename or a list of 2 fastq filenames
        :param database: the path to a valid Kraken database
        :param threads: number of threads to be used by Kraken
        :param output: output filename of the Krona HTML page

        :param return:

        """
        self._devtools = DevTools()
        self._devtools.check_exists(database)

        self.database = database
        self.threads = threads

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

        for this in self.fastq:
            self._devtools.check_exists(database)

    def run(self):
        self.kraken_output = TempFile()

        params = {
            "database": self.database,
            "thread": self.threads,
            "file1": self.fastq[0],
            "kraken_output": self.kraken_output.name,
            }

        if self.paired:
            params["file2"] = self.fastq[1]

        command = "kraken -db %(database)s %(file1)s "
        if self.paired:
            command += " %(file2)s --paired"
        command += " --threads %(thread)s --out %(kraken_output)s"

        command = command % params
        # Somehow there is an error using easydev.execute with pigz
        from snakemake import shell
        shell(command)


class KrakenDownload(object):
    """Utility to download Kraken DB and place them in a local directory

    ::

        from sequana import KrakenDownload
        kd = KrakenDownload()
        kd.download('toydb')
        kd.download('minikraken')

    A large database (8Gb) is available on synapse and has the following DOI::

        doi:10.7303/syn6171000

    It can be downloaded manually or if you have a Synapse login
    (https://www.synapse.org), you can use::

        from sequana import KrakenDownload
        kd = KrakenDownload()
        kd.downloaded("sequana_db1")
    """
    dv = DevTools()
    def download(self, name, verbose=True):
        if name == "minikraken":
            self._download_minikraken(verbose=verbose)
        elif name == "toydb":
            self._download_kraken_toydb(verbose=verbose)
        elif name == "sequana_db1":
            self._download_sequana_db1(verbose=verbose)
        else:
            raise ValueError("name must be toydb or minikraken, or sequana_db1")

    def _download_kraken_toydb(self, verbose=True):
        """Download the kraken DB toy example from sequana_data into
        .config/sequana directory

        """
        dv = DevTools()
        base = sequana_config_path + os.sep + "kraken_toydb"
        taxondir = base + os.sep + "taxonomy"
        dv.mkdir(base)
        dv.mkdir(taxondir)

        baseurl = "https://github.com/sequana/data/raw/master/"

        # download only if required
        if verbose:
            print("Downloadind the database")
        wget(baseurl + "kraken_toydb/database.idx", base + os.sep + "database.idx")
        wget(baseurl + "kraken_toydb/database.kdb", base + os.sep + "database.kdb")
        wget(baseurl + "kraken_toydb/taxonomy/names.dmp",
             taxondir +  os.sep + "names.dmp")
        wget(baseurl + "kraken_toydb/taxonomy/nodes.dmp",
            taxondir + os.sep + "nodes.dmp")

    def _download_minikraken(self, verbose=True):
        dv = DevTools()
        base = sequana_config_path + os.sep + ""
        taxondir = base + os.sep + "taxonomy"
        dv.mkdir(base)
        dv.mkdir(taxondir)
        if verbose:
            print("Downloading minikraken (4Gb)")
        wget("https://ccb.jhu.edu/software/kraken/dl/minikraken.tgz",
             base + os.sep + "minikraken.tgz")
        # unzipping. requires tar and gzip

    def _download_from_synapse(self, synid, target_dir):
        from synapseclient import Synapse
        try:
            self._synapse.get(synid, downloadLocation=target_dir)
        except:
            self._synapse = Synapse()
            self._synapse.login()
            self._synapse.get(synid, downloadLocation=target_dir)

    def _download_sequana_db1(self, verbose=True):
        dbname = "sequana_db1"
        from easydev import md5
        from sequana import sequana_config_path
        dir1 = sequana_config_path + os.sep + dbname
        dir2 = dir1 + os.sep + "taxonomy"
        self.dv.mkdir(dir1)
        self.dv.mkdir(dir2)

        print("Downloading about 8Gb of data (if not already downloaded) from"
            " Synapse into %s" % dir1) 

        from os.path import exists
        filename = dir1 + "ena_list.txt"
        if exists(filename) and md5(filename) == "a9cc6268f3338d1632c4712a412593f2":
            pass
        else:
            self._download_from_synapse('syn6171700', dir1)

        # database.idx
        filename = dir1 + "database.idx"
        if exists(filename) and md5(filename) == "2fa4a99a4f52f2f04c5a965adb1534ac":
            pass
        else:
            self._download_from_synapse('syn6171017', dir1)

        # database.kdb ; this one is large (8Gb)
        filename = dir1 + "database.kdb"
        if exists(filename) and md5(filename) == "ff698696bfc88fe83bc201937cd9cbdf":
            pass
        else:
            self._download_from_synapse('syn6171107', dir1)

        # Then, the taxonomy directory
        filename = dir1 + "names.dmp"
        if exists(filename) and md5(filename) == "10bc7a63c579de02112d125a51fd65d0":
            pass
        else:
            self._download_from_synapse('syn6171286', dir2)

        filename = dir1 + "nodes.dmp"
        if exists(filename) and md5(filename) == "a68af5a60434e2067c4a0a16df873980":
            pass
        else:
            self._download_from_synapse('syn6171289', dir2)

        filename = dir1 + "taxons.txt"
        if exists(filename) and md5(filename) == "e78fbb43b3b41cbf4511d6af16c0287f":
            pass
        else:
            self._download_from_synapse('syn6171290', dir2)
        print('done. You should have a kraken DB in %s' % dir1)










