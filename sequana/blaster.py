import os

from bioservices import NCBIblast, easyXML, ENA
import pandas as pd
from bs4 import BeautifulSoup
from easydev import shellcmd, TempFile, Progress

from sequana import FastQ, FastA
from easydev.profiler import do_profile
import numpy as np

import pysam

import json
import re
import math

# http://stackoverflow.com/questions/8730119/retrieving-json-objects-from-a-text-file-using-python
FLAGS = re.VERBOSE | re.MULTILINE | re.DOTALL
WHITESPACE = re.compile(r'[ \t\n\r]*', FLAGS)
class ConcatJSONDecoder(json.JSONDecoder):
    """Ease the reading of json-like file returned by blast when several
sequences are provided. Indeed, in such case the output is not a valid json bu a
concatenation of json

    Usage ::

    json.loads(res, cls=ConcatJSONDecoder)

    """
    def decode(self, s, _w=WHITESPACE.match):
        s_len = len(s)

        objs = []
        end = 0
        while end != s_len:
            obj, end = self.raw_decode(s, idx=_w(s, end).end())
            end = _w(s, end).end()
            objs.append(obj)
        return objs


class XML(object):
    def __init__(self, xml):
        if isinstance(xml, str):
            try: # if xml is a filename, we open
                xml = open(xml, 'r').read()
            except: # an XML string already (or easyXML instance)
                pass
            self.xml = BeautifulSoup(xml, "lxml")
        else:
            # sub strings of hits
            self.xml = "\n".join([str(x) for x in xml])


class BlastAlignment(dict):
    """A data structure to hold results from Blast alignment"""
    def __init__(self, record):
        self.record = record
        # This class should work all the time but we need to agree on one
        # single set of keys that is the ones used by EBI blast

        items = record.findAll()

        # Blast standalone 2.2.31
        if items[0].name.startswith('hsp_'):
            ncbi_mapping = {
                 #'hsp_align-len': '101',  # not stored in EBI XML
                 'hsp_bit-score': 'bits',
                 'hsp_score': 'score',
                 'hsp_evalue': 'expectation',
                 'hsp_gaps': 'gaps',
                 'hsp_hit-frame': 'frame',
                 'hsp_hit-from': 'matchseq_end',
                 'hsp_hit-to': 'matchseq_start',
                 'hsp_hseq': "matchseq",
                 'hsp_identity': 'identity',
                 'hsp_midline': 'pattern',
                 #'hsp_num': '1',  # probably the alignment attribute 'number' in EBI XML
                 'hsp_positive': 'positive',
                 'hsp_qseq': "queryseq",
                  #'hsp_query-frame': '1',  # cannot find equivalent in EBI XML
                 'hsp_query-from': 'queryseq_start',
                 'hsp_query-to': 'queryseq_end'}
            for item in items:
                try: self[ncbi_mapping[item.name]] = item.text
                except: pass
                self['strand'] = "unknown"
        else:
            # EBI case
            for item in items:
                if item.name in ['queryseq', 'matchseq']:
                    self[item.name+"_start"] = item.attrs['start']
                    self[item.name+"_end"] = item.attrs['end']
                    self[item.name] = item.text
                else:
                    self[item.name] = item.text

    def __str__(self):
        metadata = self.copy()
        pattern = self['pattern']
        metadata['ident1'] = pattern.count("|")
        metadata['ident2'] = len(pattern)
        metadata['ident_pc'] = 100.* pattern.count("|") / float(len(pattern))
        metadata['ident_pc'] = int(metadata["ident_pc"]*100)/100.
        metadata['gap1'] = "X"
        metadata['gap2'] = len(pattern)
        metadata['gap_pc'] = "X"
        metadata['q1'] = self['queryseq_start']
        metadata['q2'] = self['queryseq_end']
        metadata['m1'] = self['matchseq_start']
        metadata['m2'] = self['matchseq_end']

        txt = """
Score = %(bits)s bits (%(score)s),  Expect = %(expectation)s
Identities = %(ident1)s/%(ident2)s (%(ident_pc)s%%), Gaps = %(gap1)s/%(gap2)s (%(gap_pc)s%%)
Strand=%(strand)s

Query\t%(q1)s\t%(queryseq)s\t%(q2)s
\t\t%(pattern)s
Sbjct\t%(m1)s\t%(matchseq)s\t%(m2)s
"""  % metadata
        return txt



class BlastHit(dict):
    """A data structure to hold results of a blast hit"""
    def __init__(self, record, sequence_name="undefined"):
        self.record = record
        self.name = sequence_name
        try:
            self._ebi_case()
        except:
            self._blast_case()

    def _blast_case(self): # DOCTYPE BlastOutput
        description = self.record.find('hit_def').text
        lhs, description = description.split(" ", 1)
        db, id_, ac = lhs.split("|")
        self['database'] = db
        self['ac'] = ac
        self['description'] = description
        self['id'] = id_
        self['length'] = self.record.find('hit_len').text
        self['number'] = self.record.find('hit_num').text
        self['alignment'] = BlastAlignment(self.record.find('hsp'))
        self['sequence_name'] = self.name

    def _ebi_case(self):  # EBIApplicationResult from NCBIblast
        self['ac'] = self.record.attrs['ac']
        self['database'] = self.record.attrs['database']
        self['description'] = self.record.attrs['description']
        self['id'] = self.record.attrs['id']
        self['length'] = self.record.attrs['length']
        self['number'] = self.record.attrs['number']
        self['alignment'] = BlastAlignment(self.record.find('alignment'))
        self['sequence_name'] = self.name

    def __str__(self):
        # TODO
        txt = """
>lcl|%(database)s:%(id)s %(ac)s %(description)s Length=%(length)s
""" % self.record.attrs
        txt += self['alignment'].__str__()
        return txt


class BlastHits(XML):
    """Get Hits from a Blast XML search

    Works with a local blast, a result from NCBIblast (EBI) either
    from a filename or a string.

    ::

        xml = open("test_blast_ebi.xml").read()
        hits = BlastHits(xml)
        # equivalent to
        hits = BlastHits('test_blast_ebi.xml')

    """
    def __init__(self, xml):
        """

        xml is either a filename or xml

        """
        if isinstance(xml, str):
            if os.path.exists(xml):
                xml = open(xml, 'r').read()
            self.xml = BeautifulSoup(xml, "lxml")
        else:
            # sub strings of hits
            self.xml = "\n".join([str(x) for x in xml])

    def _get_hits(self):
        try:
            # XML input
            hits = [BlastHit(this) for this in self.xml.findAll('hit')]
        except:
            # sub-string of hits
            hits = [BlastHit(BeautifulSoup(str(hit))) for hit in self.xml]
        return hits
    hits = property(_get_hits)


class BlastResults(object):
    """Identify Hits from XML search

    Search from EBI has only 1 sequence (submission restriction)
    """
    #@do_profile()
    def __init__(self, xml, sequence_name="undefined"):

        # Figure out whether it is EBI (one query) or NCBI (several queries)
        if isinstance(xml, str):
            if os.path.exists(xml):
                xml = open(xml, 'r').read()
            self.xml = BeautifulSoup(xml, "xml")

        if self.xml.find('Iteration'):
            self.mode = "ncbi"
            self.iterations = self.xml.findAll('Iteration')
        else:
            self.mode = "ebi"

        #self.iterations = self.xml.findAll('iteration')
        self.hits = BlastHits(xml).hits
        if len(self.hits) == 0:
            columns = ['ac', 'bits', 'database', 'description',
                'expectation', 'frame', 'gaps', 'id', 'identity',
                'length', 'matchseq', 'matchseq_end', 'matchseq_start',
                'number', 'pattern', 'positives', 'queryseq',
                'queryseq_end', 'queryseq_start', 'score', 'sequence_name',
                'strand']
            self.df = pd.DataFrame(columns=columns)

        if self.mode == "ncbi":
            self.names = []
            for i, iteration in enumerate(self.iterations):
                name = iteration.findAll('Iteration_query-def')[0].text
                N = len(iteration.findAll("Hit"))
                self.names.extend([name]*N)

        self.data = []
        for i, hit in enumerate(self.hits):
            entry = hit['alignment'].copy()
            for key in ['ac', 'database', 'description', 'id', 'length',
                'number', "sequence_name"]:
                entry[key] = hit[key]
                if self.mode == "ebi":
                    entry["sequence_name"] = sequence_name
                else:
                    entry["sequence_name"] = self.names[i]
            self.data.append(entry)
        self.df = pd.DataFrame(self.data, dtype=float)

    def _get_best(self):
        if self.mode == "ebi":
            return self.df.ix[0]
        else:
            # sort by groups of queries
            pass
    best = property(_get_best)


class BlatResults(object):
    def __init__(self, blatdata):
        self.data = blatdata

    def parse(self):
        pass


class BlastResultsJson(object):

    def __init__(self, blastdata):
        self.data = blastdata
        self.json = json.loads(blastdata, cls=ConcatJSONDecoder)

    def parse(self):
        self.datalist = []
        for name, obj in enumerate(self.json):
            hits = obj['BlastOutput2']['report']['results']['search']['hits']

            for i, hit in enumerate(hits):
                hsp = hit['hsps'][0] # do we have more than one if so what do we do ?
                """if len(hit['hsps']) >1:
                    print("found %s hsp in read %s" % (len(hit['hsps']),    name   ))
                if len(hit['hsps'])>100:
                    print("  -- " + ac)
                """
                search  = obj['BlastOutput2']['report']['results']['search']
                description = hit['description'][0]
                length = hit['len']
                number = hit['num']
                lhs, description = description['title'].split(" ", 1)
                try:
                    db, id_, ac = lhs.split("|")
                except:
                    db, id_, ref, ac, _ = lhs.split("|")

                entry = {}
                #entry['nhsps'] = len(hit['hsps'])
                entry["identifier"] = lhs
                entry['database'] = db
                entry['ac'] = ac
                entry['description'] = description
                entry['id'] = id_
                #entry['length'] = length
                #entry['number'] = number
                #entry['sequence_name'] = None
                entry['bits'] = hsp['bit_score']
                entry['expectation'] = hsp['evalue']
                #entry['score'] = hsp['score']
                entry['gaps'] = hsp['gaps']
                #entry['hit_frame'] = None
                #entry['matchseq_start'] = hsp['hit_from']
                #entry['matchseq_end'] = hsp['hit_to']
                #entry["pattern"] = hsp['midline']
                entry['identity']  = hsp['identity']
                #entry['queryseq'] = hsp['qseq']
                #entry['positive'] = None
                entry['queryseq_start'] = hsp['query_from']
                entry['queryseq_end'] = hsp['query_to']
                #entry['strand'] = hsp['query_strand']
                #entry['hit_strand'] = hsp['hit_strand']
                #entry['number'] = hsp['num']
                entry['sequence_name'] = "%s" % name
                entry['taxid'] = None
                self.datalist.append(entry)

    

        self.df = pd.DataFrame(self.datalist, dtype=float)




class BlastRunOnline(object):
    """

    .. plot::

        from sequana.blaster import BlastRunOnline

        bro = BlastRunOnline()
        results = bro.run("CTTACCTTCGCATCAAGAGG")

        results.df.score.plot(marker="o")
        from pylab import xlabel, ylabel, grid
        grid(True)
        xlabel('hits', fontsize=16)
        ylabel('score', fontsize=16)


    """
    def __init__(self):
        self.kargs = {}
        self.ncbi = NCBIblast(verbose="ERROR")

    def run(self, sequence, database="em_rel_vrl"):
        # Sometines it fails, need to rerun it but we should check
        jobid = 400
        trial = 0
        while jobid == 400 and trial < 10:
            jobid = self.ncbi.run(program="blastn",
                sequence=sequence,
                stype="dna", database=database,
                email="thomas.cokelaer@pasteur.fr", **self.kargs)
        if jobid == 400:
            return []
        self.ncbi.wait(jobid)
        self.results = self.ncbi.getResult(jobid, "xml")
        # The returned object is eqsyXML from bioservices, but let us get the
        # xml data instead
        self.blast_results = BlastResults(self.results.data)
        return self.blast_results


def __runblast(sequence, database="em_rel_vrl"):
    import time
    import random
    # avoids to launch all jobs at the same time
    time.sleep(0.2*random.random())
    blaster = BlastRunOnline()
    jobid = blaster.run(sequence, database)
    return (sequence, jobid.df)


def __runblat(sequence, database, tag):
    import time
    import random
    # avoids to launch all jobs at the same time
    time.sleep(0.2*random.random())
    blaster = LocalBlat(database, tag)
    df = blaster.run(sequence)
    return (df)


def blast_multireads_online(sequences, database="em_rel_vrl", n_jobs=5):
    """This function can run blast in parallel using NCBIBlast

    !! not for production but testing

    ::

        from sequana import blaster
        from sequana import fastq
        f = fastq.FastQ("Hm2_GTGAAA_L005_R1_001.fastq.gz")
        sequences = [next(f)['sequence'].decode("utf-8") for i in range(10)]
        bl = blaster.blast_multireads(sequences)

    with 10 jobs in // takes about 200 seconds for 200 sequences on the
    em_rel_vrl database

    :return:  list of tuples. First element is the sequence and second item is
    the dataframe from BlastResults
    """
    from easydev import MultiProcessing
    mc = MultiProcessing(n_jobs)
    for sequence in sequences:
        mc.add_job(__runblast, sequence, database)
    mc.run()
    return mc.results


def blat_multireads(sequences, database, njobs=4):
    from easydev import MultiProcessing
    mc = MultiProcessing(njobs)
    for i,sequence in enumerate(sequences):
        mc.add_job(__runblat, sequence, database,tag=i+1)
    mc.run()
    df = pd.concat(mc.results)
    return df


class LocalBlat(object):

    """

    PSL format expained here: http://www.ensembl.org/info/website/upload/psl.html#tracklines

    Fields are space-separated, and all 21 are required.

    matches - Number of matching bases that aren't repeats.
    misMatches - Number of bases that don't match.
    repMatches - Number of matching bases that are part of repeats.
    nCount - Number of 'N' bases.
    qNumInsert - Number of inserts in query.
    qBaseInsert - Number of bases inserted into query.
    tNumInsert - Number of inserts in target.
    tBaseInsert - Number of bases inserted into target.
    strand - defined as + (forward) or - (reverse) for query strand. 
            In mouse, a second '+' or '-' indecates genomic strand.
    qName - Query sequence name.
    qSize - Query sequence size.
    qStart - Alignment start position in query.
    qEnd - Alignment end position in query.
    tName - Target sequence name.
    tSize - Target sequence size.
    tStart - Alignment start position in query.
    tEnd - Alignment end position in query.
    blockCount - Number of blocks in the alignment.
    blockSizes - Comma-separated list of sizes of each block.
    qStarts - Comma-separated list of start position of each block in query.
    tStarts - Comma-separated list of start position of each block in target.

    Blast 
    Blast 8 and 9 have those fields
    # Fields: Query id, Subject id, % identity, alignment length, mismatches,
    # gap openings, q. start, q. end, s. start, s. end, e-value, bit score

    Blast9as comapred to blast8 separated hits for each entry with the query
    name (redudnant), database and header

    """

    def __init__(self, dbname, tag=1):
        self.dbname = dbname
        self.params = {
            "output": "output_%s.psl" % tag,
            "outfmt": "blast8",
            "database": self.dbname}
        # blast unlinke blast does not give back the description...
        # so we need to get it back ourself from the database if provided
        self.fastaDB = FastaDB(self.dbname)


    def run(self, sequence):
        cmd = "blat %(database)s %(sequence)s %(output)s -out=%(outfmt)s -noHead"

        fh = TempFile()
        ftemp = open(fh.name, "w")

        if isinstance(sequence, str):
            sequence = [sequence]

        for i, this in enumerate(sequence):
            ftemp.write(">{}\n".format(i+1))
            ftemp.write(this+"\n")
        ftemp.close()
        self.params['sequence'] = fh.name
        cmd = cmd % self.params
        shellcmd(cmd)
        fh.delete()

        # Fields: Query id, Subject id, % identity, alignment length, mismatches, gap
        # openings, q. start, q. end, s. start, s. end, e-value, bit score

        # We change the naming to agree with those from blast
        columns = ["query_id", "identifier", "identity", "alignment_length", 
            "mismatches", "gaps", "queryseq_start", "queryseq_end","s_start",
            "s_end", "expectation","bits"]   

        try:
            df = pd.read_csv(self.params['output'], sep='\t', header=None)
        except:
            # It is possible that the output is empty
            df = pd.DataFrame()
        df.columns = columns

        def _get_acid(lhs):
            try:
                db, id_, ac = lhs.split("|")
            except:
                db, id_, ref, ac, _ = lhs.split("|")
            return db, id_, ac
        acs = [_get_acid(this)[2] for this in df.identifier]
        ids = [_get_acid(this)[1] for this in df.identifier]
        dbs = [_get_acid(this)[0] for this in df.identifier]
        df["id"] = ids
        df["ac"] = acs
        df["database"] = dbs

        # now, each query may have several hits. We want to keep only one row at
        # max per each query. In Pandas, we group by bits before hand so that we 
        # can then later the first element in each group
        # In principle the output in blat is already sorted but we sort again 
        df = df.sort_values('bits', ascending=False).groupby('query_id', as_index=False)
        df = df.first()
        # query id was just from blat and is not required anymore

        df.set_index("query_id",  inplace=True)
        df.drop("s_start", axis=1, inplace=True)
        df.drop("s_end", axis=1, inplace=True)

        descriptions = []
        for this in df.identifier.values:
            de = self.fastaDB.df[self.fastaDB.df.name == this].comment.iloc[0]
            descriptions.append(de)
        df['description'] = descriptions

        # Finally, some sequences may not have been found..
        for i, this in enumerate(sequence):
            if i+1 not in df.index:  #Note that query_id starts at 1 hence +1
                df.ix[i+1] = [None] * len(df.columns)

        df = df.sort_index().reset_index(drop=True)

        return df


class LocalBlast(object):
    def __init__(self, dbname, tag=1):
        self.dbname = dbname
        self.tag = str(tag)

        self.params = {
            "database": self.dbname,
            "outfmt": 13,
            "expectation":10,
            "max_hsps": 10,
            "max_target_seqs":1,
            "output": "out_%s.xml" % tag}

        #outfmt : 5 for xml, 13 for json

    def run(self, sequence):
        if self.tag:
            output = "out_%s.xml" % self.tag
        # -dust no and -task blastn are used to mimic results from EBI
        cmd = "blastn -task blastn -num_threads 8 -query %(sequence)s  -db %(database)s "
        cmd += " -evalue %(expectation)s "
        cmd += " -max_hsps %(max_hsps)s "
        cmd += " -max_target_seqs %(max_target_seqs)s -dust no"
        cmd += " -outfmt %(outfmt)s > %(output)s "
        """
        -num_alignments 50 -reward 1 -penalty -3 -gapopen 5 -gapextend 2

        """
        fh = TempFile()
        ftemp = open(fh.name, "w")
        if isinstance(sequence, str):
            sequence = [sequence]
        for i, this in enumerate(sequence):
            ftemp.write(">{}\n".format(i+1))
            ftemp.write(this+"\n")
        ftemp.close()
        self.params['sequence'] = fh.name
        cmd = cmd % self.params
        shellcmd(cmd)
        fh.delete()
        with open(self.params["output"]) as fin:
            results = fin.read()

        bjson = BlastResultsJson(results)
        bjson.parse()
        df = bjson.df
        # Finally, some sequences may not have been found..

        df.set_index("sequence_name", inplace=True)
        for i, this in enumerate(sequence):
            if i not in df.index:  #Note that query_id starts at 1 hence +1
                df.ix[i] = [None] * len(df.columns)
        df = df.sort_index().reset_index(drop=True)
        return df


class FastaDB(object):
    """Handler for a FASTA file

    Can add a new fasta sequence. Checks if it is present otherwise adds it.
    Can check the names of the identifier, or sequence length within a
    dataframe. Can also create a blast database.


    """
    def __init__(self, filename="db/measles.fasta"):
        self.filename = filename
        if os.path.exists(self.filename) is False:
            with open(self.filename, "a") as fin:
                pass
                # creates empty file
        self._init()

    def _init(self):
        # Let us store the name and description and length of the sequence
        f = FastA(self.filename)

        from easydev import Progress
        pb = Progress(len(f))
        data = []
        for i, this in enumerate(f):
            data.append((this.name, this.comment, len(this.sequence)))
            pb.animate(i+1)
        print()
        self.df = pd.DataFrame.from_records(data,
            columns=['name', 'comment', 'seq_length'])

        def _getids(x):
            try:
                return x.split("|")[1]
            except:
                return x
        self.df['id'] = self.df['name'].apply(lambda x: _getids(x))

    def append_fasta(self, fasta):
        # check that this fasta is not already in the DB
        # If not, add to the local dataframe
        header, sequence = fasta.split("\n", 1)
        name, comment = header.split(" ")
        if name.startswith('>'):
            name = name[1:]
        if "|" in name:
            id_ = name.split("|")[1]
        else:
            id_ = name
        length = len(sequence.replace('\n', ''))

        if id_ in self.df['id'].values:
            print('%s already present' % id_)
        else:
            with open(self.filename, "a") as fh:
                fh.write(fasta)
                if fasta.endswith('\n') is False:
                    fh.write("\n")
            self.df.append(
                {   "name":name,
                    "comment":comment,
                    "seq_length":length,
                    "id":id_}, ignore_index=True)

    def updateDB(self):
        # update the local DB for blast
        cmd = "makeblastdb -in %s -dbtype nucl" % self.filename
        shellcmd(cmd)

    def get_stats(self):
        stats = {}
        stats['nseq'] = len(self)
        f = pysam.FastxFile(self.filename)
        stats['nbp'] = sum([len(this.sequence) for this in f])
        return stats

    def __str__(self):
        txt = "FastaDB with %s sequences" % len(self)
        return txt

    def __len__(self):
        return len(self.df)


class Contaminant(object):
    """

    Input file name must be a valid FastQ



    http://www.ebi.ac.uk/Tools/sss/ncbiblast/help/index-nucleotide.html

    EST: expressed sequence tag
    GSS: Genome survey sequence
    HTC: High throughput cDNA
    HTG: High throughput genomic
    WGS/ whole genome shotgun
    STS: 
    TSA: transcriptomic shotgun assembly

    em_rel: all
    em_rel_env: environment
    em_rel_fun: fungi
    em_rel_hum: human
    em_rel_inv: invertebrate
    em_rel_mam: mammal without humans and rodents
    em_rel_mus: mouse mus musculus
    em_rel_phg: phages
    em_rel_pln: plants
    em_rel_rod: rodents but not mouse
    em_rel_pro: prokaryote
    em_rel_syn: synthetic constructs
    em_rel_tgn: transgenic constructs
    em_rel_unc: unspecified
    em_rel_vrl: viruses
    em_rel_vrt: vrtebrates excluding human, mouse, rodents

    """
    def __init__(self, fastq_filename, database="em_rel",
            local_fasta_db="db/measles.fasta", online=True, fastaDB=False):
        self.ena = ENA()

        self.chunk = 10  # number of sequences in a single chunk to look at
        self.online = online

        # allow to run blast locally
        self.local_blast = LocalBlast(local_fasta_db)
        self.local_blat = LocalBlat(local_fasta_db)
        self.alignment = "blat"

        # The fastq to process
        if "fastq" in fastq_filename:
            self.fastq = FastQ(fastq_filename)
        elif "fasta" in fastq_filename:
            self.fastq = FastA(fastq_filename)

        self.database = database
        self.local_fasta_db = local_fasta_db

        if fastaDB:
            self.fastaDB = FastaDB(filename=local_fasta_db)

        self._indices = []
        self.sequences = []
        self.identifiers = []
        self.search_missing = False
        self._nreads = 0
        # rough approximation of the threshold to consider good and bad hits
        # the siwe of the query may changed, let us be conservative and pick up
        # the largest one
        try:
            n = max([len(x['sequence']) for x in self.fastq])
        except:
            n = max([len(x.sequence) for x in self.fastq])

        if fastaDB:
            m = self.fastaDB.get_stats()['nbp']
            self.bits_threshold = math.ceil(math.log2(m*n*len(self.fastq)))

            self.verbose = True
            print("Threshold for bits score set to %s" % self.bits_threshold)
        else:
            self.bits_threshold = 50

    def run(self, chunk=2000):
        # 350 seconds using the virus DB on 110000 reads settings search_missing
        # to False
        import time
        t1 = time.time()
        self.chunk = chunk
        N = int(len(self.fastq)/float(chunk)) + 1
        pb = Progress(N)
        for i in range(N):
            self.analysis()
            pb.animate(i+1)
        t2 = time.time()
        if self.verbose:
            print(t2-t1)

    def run_blat(self, chunk=1000, njobs=4, sequences=None):
        """

        Takes chunk of reads (e.g., 1000), uses the blat_multiread function
        ito use N threads (njobs). blast_multireads  uses local database
        reference to run blat on it.

        The local DB should be handled by users or build using 
        run_blast, which uses an online search on missing reads.


        """
        sequences = [this['sequence'].decode("utf-8") for this in self.fastq]
        from easydev import split_into_chunks
        if len(sequences) < 1000:
            nchunks = 4
        else:
            nchunks = int(len(sequences)/1000.) + 1
        this = split_into_chunks(sequences, nchunks)
        res = blat_multireads(list(this), database=self.local_fasta_db, njobs=njobs)
        res.reset_index(drop=True, inplace=True)
        self.results = res
        return res

    def get_taxon_ena(self, identifier):
        data = self.ena.get_data(identifier, frmt="text")
        taxon = [x for x in data.decode().split("\n") if "taxon" in x]
        taxon = taxon[0].split('taxon:')[1]
        taxon = taxon.replace("\"", "")
        return taxon

    def get_taxons_eutils(self):
        # fill the taxon ids
        # FIXME. use accession for now but need to check if we do not want the
        # identifiers instead ? WE DO indeed it seems the id links to a taxon
        # all the time but not necceseraly the accession. e.g 
        df = self.df_results_with_dummies()

        # There are 2 types of identifiers. Those from GI are pure identifiers
        # with numbers; 
        ids =  list(set(df.id.dropna()))
        gi_ids = []
        ena_ids = []
        for this in ids:
            try:
                int(this)
                gi_ids.append(int(this))
            except:
                ena_ids.append(this)
        # First, we deal with the GI ids
        taxons = get_taxon_from_ids(gi_ids)
        taxon_dict = dict([(k,v) for k,v in zip(gi_ids, taxons)])
        # then the ENA ones
        for this in ena_ids:
            taxon = self.get_taxon_ena(this)
            taxon_dict[this] = taxon
        return taxon_dict

    def analysis(self):
        reads = self.get_chunk_of_reads()
        sequences = [r['sequence'].decode('utf-8') for r in reads]
        identifiers = [r['identifier'].decode('utf-8') for r in reads]
        if len(sequences) == 0:
            if self.verbose:
                print('no more sequences to process')
            return

        # Here, the N sequences are named internally 1,2,...N
        if self.alignment == "blast":
            df = self.local_blast.run(sequences)
        else:
            df = self.local_blat.run(sequences)

        # We may want to cross out some values
        assert len(df) == len(sequences) 
        badhits = sum(df.bits < self.bits_threshold)
        df[df.bits < self.bits_threshold] = None
        df.fillna(np.nan, inplace=True)

        try:
            self.results = self.results.append(df)
            self.results.reset_index(inplace=True, drop=True)
        except:
            self.results = df

    def filter_out_low_quality(self):
        badhits = sum(self.results.bits < self.bits_threshold)
        self.results[self.results.bits < self.bits_threshold] = None
        self.results.fillna(np.nan, inplace=True)

    def update_missing(self, N=None):
        if self.verbose:
            print("Running online blast ")

        self._indices = self.results.index[self.results.bits.isnull()]
        self._sequences =  [x['sequence'].decode() for i,x in enumerate(self.fastq) 
                            if i in self._indices]
        if N is not None:
            self._indices = self._indices[0:N]
            self._sequences = self._sequences[0:N]

        self.res = blast_multireads_online(self._sequences, database=self.database)
        if self.verbose:
            print('Fetching new FASTA')
        for i, this in enumerate(self.res):
            seq, df = this
            if len(df) ==0 :
                print("No hits for %s... in %s . Use another DB" % (seq[0:20], self.database))
            else:
                self.results.ix[self._indices[i]] = df.ix[0]
                # found it so, let us now download the fasta
                # and update the local DB
                identifier = df.ix[0]['id']
                #self.add_fasta_in_db(identifier)

    def add_fasta_in_db(self, identifier):
        if identifier not in self.fastaDB.ids:
            if self.verbose:
                print("Downloading %s from ENA and updating local DB" % identifier)
            fasta = self.download_fasta(identifier)
            self.fastaDB.append_fasta(fasta)
        else:
            if self.verbose:
                print("%s already in DB " % identifier)
        self.fastaDB.updateDB()

    def get_chunk_of_reads(self):
        # must use next not a loop since the loop restarts at position 0
        # but we also must keep track of the position otherwise we will loop
        # again...
        reads = []
        count = 0
        N = self.chunk
        while count < self.chunk and self._nreads <len(self.fastq):
            try:
                reads.append(next(self.fastq))
            except:
                pass # we reached the end of the sequences
            count += 1
            self._nreads +=1
        return reads

    def get_lineage(self, identifier):
        """
        :param identifier: a valid ENA accession
        """
        data = self.ena.get_data(identifier, "txt")
        lineage = ""
        for line in data.decode().split("\n"):
            if line.startswith("OS"):
                name = line.split(" ", 1)[1].strip()
            if line.startswith('OC'):
                lineage += line.split(" ", 1)[1].strip()
        return name, lineage

    def get_lineage_from_taxons(self, taxons):
        from bioservices import EUtils
        e = EUtils()
        res = e.EFetch("taxonomy", taxons)
        self.xml = easyXML(res)

        # There is one lineage per taxon so we can use the findAll
        lineage = [x.text for x in self.xml.findAll("lineage")]

        # We now want to get the scientific name for each taxon. Note,ever
        # that there are several taxon children-tags within a taxon each hving
        # a scientific name, so we first use getchildren() to make sure toloop
        # over a the children only
        scnames = []
        for taxon in self.xml.root.getchildren():
            name = taxon.findall("ScientificName")[0].text
            scnames.append(name)

        return lineage, scnames

    def download_fasta(self, name):
        fasta = self.ena.get_data(name, "fasta")
        fasta = fasta.decode('utf-8')
        return fasta

    def plot(self):
        res = self.results.groupby("description").size()
        res = res.reset_index().sort_values(by=0)
        res.columns = ["description", "count"]
        return res

    def __str__(self):
        txt = """
        Number of input reads: %s
        Number of analysed reads: %s
        Number of found reads %s
        Number of unclassified hits reads %s
        """
        N = len(self.fastq)
        NA = len(self.results)
        C = len(self.results)
        U =  NA-C
        return txt % (N, NA, C, U)

    def save_classified(self, tag):
        """
        # SAme in same kinf of format as kraken:
        C identifier1 taxid length quality
        U identifier2 0 length quality

        C= classified
        U = unclassified
        taxid to be provided if classified, otherwise 0

        """
        filename_class = "classified_%s.fastq" % tag
        filename_unclass = "unclassified_%s.fastq" %tag
        with open(filename_class, "wb") as foutc:
            with open(filename_unclass, "wb") as foutu:

                for i, reads in enumerate(self.fastq):
                    if isinstance(self.results[i], str):
                        foutu.write(reads['identifier']+b"\n")
                        foutu.write(reads['sequence']+b"\n")
                    else:
                        foutc.write(reads['identifier']+b"\n")
                        foutc.write(reads['sequence']+b"\n")

    def build_krona(self, output_filename):
        # Some accession are ENA and others from NCBI GI ID

        df = self.df_results_with_dummies()
        counts = df.groupby('id').size()
        #counts.drop(0, inplace=True)
        counts.sort_values(ascending=False, inplace=True)
        counts = counts.to_frame()
        counts.columns=['count']

        taxons_dict = self.get_taxons_eutils()
        taxons = list(taxons_dict.values())
        counts['taxid'] = [taxons_dict[k] for k in counts.index]


        lineage = []
        scnames = []
        tolookat = counts['taxid'].values
        for i in range(int(len(tolookat)/150.)+1):
            thisvalues = tolookat[i*150:(i+1)*150]
            lins, names = self.get_lineage_from_taxons(",".join(thisvalues))
            lineage.extend(lins)
            scnames.extend(names)

        assert len(lineage) == len(taxons)
        assert len(scnames) == len(taxons)

        N = len(self.fastq)
        NA = len(self.results)
        C = len(self.df_results)
        U =  NA-C

        self.lineage = lineage
        self.counts = counts

        #counts.index = [int(x) for x in counts.index]

        with open(output_filename, "w") as fout:
            for i, this in enumerate(lineage):
                line = str(counts['count'].iloc[i])+"\t"+"\t".join(this.split(';'))
                line += " " + scnames[i]
                fout.write(line+'\n')
            fout.write("%s\t%s" % (U, "Unclassified"))


def get_taxon_from_ids(idlist):

    from bioservices import EUtils
    eutils = EUtils(cache=True)
    pb = Progress(len(idlist))
    taxons = []
    for i, id_ in enumerate(idlist):
        taxons.append(get_taxon_from_one_id(id_, eutils))
        pb.animate(i+1)
    return taxons


def get_taxon_from_one_id(id_, eutils=None):
    """

    """
    if eutils is None:
        from bioservices import EUtils
        eutils = EUtils()
    ret = eutils.ELink("taxonomy", "nucleotide", [id_])
    parsed = eutils.parse_xml(ret, 'EUtilsParser')

    assert str(id_) == parsed.eLinkResult.LinkSet.IdList.Id
    return parsed.eLinkResult.LinkSet.LinkSetDb.Link.Id




def _get_taxon_from_gi_ids(idlist):
    """
    from bioservices import EUtils
    e = EUtils()
    ret = e.ELink("taxonomy", "nucleotide", idlist)

    # This API may change
    parsed = e.parse_xml(ret, 'EUtilsParser')

    # Note that may be not all ids are correct so we need to check that as well
    identifiers = parsed.eLinkResult.LinkSet.IdList.Id
    try:
        taxons = [x.Id for x in parsed.eLinkResult.LinkSet.LinkSetDb.Link]
    except:
        taxons = [parsed.eLinkResult.LinkSet.LinkSetDb.Link.Id]

    print(len(identifiers), len(taxons))
    results = []
    for identifier in idlist:
        try:
            index = idlist.index(identifier)
            results.append(taxons[index])
        except:
            results.append(None)
    return results
    """





