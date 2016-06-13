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
from easydev import execute, TempFile, Progress
import pylab


__all__ = ['KrakenResults', "KrakenBuilder", 
    "KrakenPipeline", "KrakenAnalysis"]



class KrakenResults(object):
    """Translate Kraken results into a Krona-compatible file


    If you run a kraken analysis, you will end up with a file e.g. kraken.out

    You could use kraken-translate but then you need extra parsing to convert
    into a Krona-compatible file. Here, we take the output from kraken and
    directly transform to the krona-compatible file.

    ::

        k = KrakenResults()
        k.kraken_to_krona()

    Then in a shell, ::


    """

    def __init__(self, filename="kraken.out", verbose=True):
        self.filename = filename
        self.verbose = verbose
        from biokit import Taxonomy
        self.tax = Taxonomy(verbose=self.verbose)
        self._data_created = False

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
            fout.write("%s\t%s" % (self.unclassified, "Unclassified"))
        self._data_created = True

    def plot(self, kind="pie", cmap="copper", threshold=1, radius=0.7, **kargs):
        """A simple non-interactive plot of taxons

        A Krona Javascript output is also available in :meth:`kraken_to_krona`
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
        from easydev import execute
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


        ka = KrakenAnalysis("R1.fastq.gz", database="krakendb")

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


