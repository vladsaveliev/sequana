

"""
IDEA

zcat test.fastq.gz | python filter ^ACGT | gzip > fitlered.fastq.gz

import re
query=re.compile("^ACGT")
if query.search(sequence):
    write(sequence)

Interesting as well for creating a generic zipfile opener

http://www.genomearchitecture.com/2014/01/how-to-gunzip-on-the-fly-with-python

"""
import time
import zlib
from itertools import islice
from easydev import do_profile
import gzip


class Identifier(object):
    def __init__(self, identifier, version=None):
        self.identifier = identifier
        if version is None:
            self.version = self._infer_version()
        else:
            self.version = version
        self.info = self._interpret()

    def _interpret(self):
        if self.version == "Illumina_1.8+":
            return self._interpret_illumina_1_8()
        elif self.version == "Illumina_1.4+":
            return self._interpret_illumina_1_8()
        else:
            return {'identifier': self.identifier[:]}

    def _infer_version(self):
        """

    old illumina # is followed by 0 for no indexing but appear to use #NNNNNN for multiplex ID::

            @HWUSI-EAS100R:6:73:941:1973#0/1

        New Illumina (e.g., 1.8)::

            @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG


        NCBI ::

            @FSRRS4401BE7HA [length=395] [gc=36.46] [flows=800] [phred_min=0] \
                                           [phred_max=40] [trimmed_length=95]

http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm

        """
        if self.identifier.count(':') == 4:
            left, right = self.identifier.rsplit(":", 1)
            if right.count("#") == 1 and right.count("\\") == 1:
                return "Illumina_1.4+"


        if self.identifier.count(" ") == 1:
            left, right = self.identifier.split(" ")
            if left.count(':') == 6 and right.count(':') == 3:
                return "Illumina_1.8+"

        return "Custom or Unknown identifier"

    def _interpret_illumina_1_8(self):
        """

        @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

        Note the space and : separators
        """
        assert self.identifier[0] == "@"
        # skip @ character
        identifier = self.identifier[1:]
        # replace spaces by : character
        identifier = ' '.join(identifier.split())
        identifier = identifier.replace(' ', ':')
        items = identifier.split(':')
        if len(items) != 11:
            raise ValueError('Niumber of items in the identifier should be 11')
        res = {}
        res['identifier'] = self.identifier[:]
        res['instrument'] = items[0]
        res['run_id'] = items[1]
        res['flowcell_id'] = items[2]
        res['flowcell_lane'] = items[3]
        res['tile_number'] = items[4]
        res['x_coordinate'] = items[5]
        res['y_coordinate'] = items[6]
        res['member_pair'] = items[7]
        res['filtered'] = items[8]
        res['self._control_bits'] = items[9]
        res['index_sequence'] = items[10]
        return res

    def _interpret_identifier_1_4(self):
        # skip @ character
        identifier = self.identifier[1:]
        identifier = identifier.replace('#', ':')
        identifier = identifier.replace('/', ':')
        items = identifier.split(':')

        if len(items) != 7:
            raise ValueError('Number of items in the identifier should be 7')
        # ['@HWUSI-EAS100R', '6', '73', '941', '1973#0/1']
        res['identifier'] = self.identifier[:]
        res['instrument_name'] = items[0]
        res['flowcell_lane'] = items[1]
        res['tile_number'] = items[2]
        res['self._x_coordinate'] = items[3]
        res['self._y_coordinate'] = items[4]
        res['self._index'] = '#' + items[5]
        res['self._member_pair'] = '/' + items[6]


    def __repr__(self):
        return "Identifier (%s)" % self.version



class FASTQ(object):
    """

    extract first 100k lines

    f = FastQ("")
    f.extract_head(100000, output='test.fastq')
    f.extract_head(100000, output='test.fastq.gz')

    equivalent to

    zcat myreads.fastq.gz | head -100000 | gzip > test100k.fastq.gz


    Counts the number of lines ::

        f.count_lines()

    or reads (assuming 4 lines per read)::

        f.count_reads()



    """
    def __init__(self, filename, verbose=False):
        self.filename = filename
        self.verbose = verbose
        # open the file in reqd mode
        self._input = open(self.filename, "r")
        self._infer_content()
        #self.count_reads()

    def _infer_content(self):
        record = self.next()
        return record

    def count_reads_gz(self, CHUNKSIZE=65536):
        # this is fast. On a 63M reads, takes 21 seconds as
        # compared to 46 s (real) and 1.13 (user) with zcat | wc
        # wc seems slow (same effects with uncompressde file).
        # Using gzip.open and reading lines is slower by a factor 10
        # recipe found http://wiki.glitchdata.com/index.php?title=Python:_File_Compression_and_Decompression
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        f = open(self.filename, 'rb')
        buf = f.read(CHUNKSIZE)
        count = 0
        while buf:
            outstr = d.decompress(buf)
            count += outstr.count("\n")
            buf = f.read(CHUNKSIZE)
        f.close()
        return count

    def count_lines(self):
        """Return number of lines


        This is 40 times faster than using SeqIO from
        """
        if self.filename.endswith("gz"):
            count = self.count_reads_gz()
            return count
        else:
            count = self._count_reads_buf()
            return count

    def count_reads(self):
        #TODO assert multiple de 4
        return self.count_lines() / 4

    """
    #slower than cound_reads_buf
    def count_read_simple(self):
        with open(self.filename, 'r') as f:
            for i, l in enumerate(f):
                pass
        i += 1
        return i
    """
    def _count_reads_buf(self, block=1024*1024):
        # 2x faster than count_reads_simple
        # 0.12 seconds to read 3.4M lines
        # surprinsingly much faster than unix command wc
        # on 2M reads, takes 0.1 seconds whereas wc takes 1.2 seconds
        f = open(self.filename, 'r')
        lines = 0
        read_f = f.read
        buf = read_f(block)
        while buf:
            # todo : use EOF universal ?
            lines += buf.count('\n')
            buf = read_f(block)
        f.close()
        return lines

    def extract_head(self, N, output_filename):
        # equivalent to
        # zcat input.fasta | head -400000 > out.fasta
        # 3 times slower than cat file | head -1000000 > test.fastq
        # but remains pretty fast (0.5 seconds for 2M reads).
        # todo: compress the file
        if self.filename.endswith(".gz"):
            self._extract_head_gz(N, output_filename)
        else:
            self._extract_head(N, output_filename)

    def _extract_head(self, N, output_filename):
        with open(self.filename, 'r') as fin:
            if output_filename.endswith("gz"):
                output_filename_nogz = output_filename.replace(".gz", "")
                with open(output_filename_nogz, 'w') as fout:
                    fout.writelines(islice(fin, N))
                print("zipping the file with gz in a shell")
                import subprocess
                s = subprocess.Popen(["gzip", output_filename_nogz, "-f"])
            else:
                with open(output_filename, 'w') as fout:
                    fout.writelines(islice(fin, N))

    def _extract_head_gz(self, N, output_filename="test.fastq.gz", level=6, CHUNKSIZE=65536):
        """

        if output not compressed, this is 20% faster than
        "zcat file | head -1000000 > output.fastq

        If output is compressed, this is equivalent to :
        "zcat file | head -1000000 | gzip > output.fastq

        Tested under Python 2.7 , Linux box.
        """
        # make sure N is integer
        N = int(N)
        # as fast as zcat file.fastq.gz | head -200000 > out.fasta
        d = zlib.decompressobj(16 + zlib.MAX_WBITS) # this is to supress the header
        f = open(self.filename, 'rb')
        buf = f.read(CHUNKSIZE)
        count = 0

        if output_filename.endswith(".gz"):
            zip_output = True
            fout = gzip.open(output_filename, "wb", compresslevel=level)
        else:
            fout = open(output_filename, "wb")
            zip_output = False

        dd = zlib.compressobj(level) # second qrg is deflated, third is zlib.MAX_WBITS

        while buf:
            outstr = d.decompress(buf)
            count += outstr.count("\n")
            if count > N:
                # there will be too many lines, we need to select a subset
                missing = count - N
                outstr = outstr.strip().split("\n")
                NN = len(outstr)
                outstr = "\n".join(outstr[0:NN-missing-1]) + "\n"
                if zip_output:
                    fout.write(outstr)
                else:
                    fout.write(outstr)
                break
            if zip_output:
                fout.write(outstr)
            else:
                fout.write(outstr)
            buf = f.read(CHUNKSIZE)

        f.close()
        fout.close()
        return count

    def extract(self, N, output_filename):
        if self.filename.endswith('gz'):
            raise NotImplementedError
        else:
            self.extract_head(N, output_filename)

    def select_random(self, N, output_filename):
        pass

    def random(self, N=10000, output_filename="test.fastq",
               bp=50, quality=40):
        """
        N here is the number of reads
        """
        # a completely random fastq
        from .phred import quality
        with open(output_filename, "w") as fh:
            count = 1
            template = "@Insilico\n"
            template += "%(sequence)\n"
            template += "+\n"
            template += "%s(quality)\n"
            fh.writelines(template % {
                'sequence': "".join(["ACGT"[random.randint(0,3)] for this in xrange(bp)]),
                'quality': "".join()})

        # quality could be q function for a distribution

    def splitting(self, N):
        pass

    def joining(self, pattern, output_filename):
        """

        zcat Block*.fastq.gz | gzip > combined.fastq.gz

        """
        pass

    def __iter__(self):
        N = self.count_reads()
        for i in xrange(0, N):
            yield self.next()

    def next(self):
        # reads 4 lines
        data = islice(self._input, 4)
        identifier = data.next()
        sequence = data.next()
        skip = data.next()
        quality = data.next()
        return {'identifier':identifier,
                'sequence': sequence,
                'quality': quality}






class FastQC(object):
    def __init__(self, filename):
        self.fastq = FASTQ(filename)

    def plot_quality(self):
        import pandas as pd
        import phred
        N = self.fastq.count_reads()
        data = []
        for i, record in enumerate(self.fastq):
            data.append(phred.Quality(record['identifier']).quality)
        return data





