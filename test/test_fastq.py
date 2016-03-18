from sequana import fastq
from . import data
pathdata = data.__path__[0]
from easydev import TempFile
import os
from nose.plugins.attrib import attr



def test_fastq_unzipped():

    # isntanciation 
    f = fastq.FastQ(pathdata + os.sep + "reads.fastq")
    assert f.data_format == "Illumina_1.8+"
    # count lines
    # rune it twice because we want to make sure re-running count_lines
    # (decompression with zlib) works when run again.
    assert f.count_lines() == 1000
    assert f.count_lines() == 1000
    assert f.count_reads() == 250
    assert f.count_reads() == 250

    # extract head of the file
    ft = TempFile() 
    f.extract_head(100, ft.name)
    fcheck = fastq.FastQ(ft.name)
    assert fcheck.count_lines() == 100
    ft.delete() 

    # extract head of the file and zip output
    ft = TempFile() 
    f.extract_head(100, ft.name)
    fcheck = fastq.FastQ(ft.name)
    assert fcheck.count_lines() == 100
    ft.delete() 


def test_fastq_zipped():

    f = fastq.FastQ(pathdata + os.sep + "reads_gz.fastq.gz")
    assert f.count_lines() == 1000
    assert f.count_reads() == 250


def test_split():
    # general tests 
    f = fastq.FastQ(pathdata + os.sep + "reads.fastq")
    try:
        f.split_lines(250) # not a multiple of 4
        assert False
    except:
        assert True

    # 
    outputs = f.split_lines(500)
    assert len(outputs) == 2 
    remove_files(outputs)
    outputs = f.split_lines(256)
    assert len(outputs) == 4
    remove_files(outputs)

    # Now tests the zip/unzip cases
    f = fastq.FastQ(pathdata + os.sep + "reads.fastq")
    outputs = f.split_lines(500, gzip=False)
    remove_files(outputs)
    outputs = f.split_lines(500, gzip=True)
    remove_files(outputs)

    f = fastq.FastQ(pathdata + os.sep + "reads_gz.fastq.gz")
    outputs = f.split_lines(500, gzip=False)
    remove_files(outputs)
    outputs = f.split_lines(500, gzip=True)
    remove_files(outputs)


def remove_files(filenames):
    for filename in filenames:
        os.remove(filename)


def test_identifiers():
    f = fastq.FastQ(pathdata + os.sep + "reads.fastq")
    identifier = fastq.Identifier(f.next()["identifier"])
    assert identifier.version == 'Illumina_1.8+'


def test_fastqc():
    qc = fastq.FastQC(pathdata + os.sep + "reads.fastq")
    qc.boxplot_quality()
    qc.histogram_gc_content()
    qc.imshow_qualities()
    qc.histogram_sequence_lengths()





