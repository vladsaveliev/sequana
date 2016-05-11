"""

Used for adapter removal benchmarks. 

May not be used on long term. 

May be removed


Lots of different formats are used to store/read adapters...

One of them, which is convenient since it exists already is Fasta::

    >NextFlex_PCR_Free_adapter1   NextFlex_PCR_Free_adapter1
    GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG

Note however, that the name is not standard here. If you want use this fasta in
other tools, it may fail, so let us add the missing bits that is 


    >NextFlex_PCR_Free_adapter1|kraken:taxid|10000001   NextFlex_PCR_Free_adapter1
    GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG

AlienTrimmer format is home-made.advantages can add comments but need another parser
cleanngs uses a two columns format for FWD and REV

"""

import pandas as pd

def fasta_fwd_rev_to_columns(file1, file2=None, output_filename=None):
    """Reads FWD and (optional) REV adapters in FASTA and save
    into a column-style file


    """

    import pysam
    f1 = pysam.FastxFile(file1)
    if output_filename is not None:
        fout = open(output_filename, "w")
    if file2:
        f2 = pysam.FastxFile(file2)
        for read1, read2 in zip(f1, f2):
            txt = "%s %s" % (read1.sequence, read2.sequence)
            if output_filename is None:
                print(txt)
            else:
                fout.write(txt+"\n")
    else:
        for read1 in f1:
            txt = "%s" % read1.sequence
            if output_filename is None:
                print(read1.sequence, read2.sequence)
            else:
                fout.write(txt+"\n")
    if output_filename is not None:
        fout.close()


def adapters_files_to_list(filename1, filename2):

    fh1 = open(filename1, 'r')
    fh2 = open(filename1, 'r')
    data1 = fh1.readlines()
    data2 = fh2.readlines()
    fh1.close()
    fh2.close()

    len(data1) == len(data2), "incompatible files. Must have same length"

    fh = open("adapters_list.fa", 'w')
    count = 0
    for line1, line2 in zip(data1, data2):
        line1 = line1.strip()
        line2 = line2.strip()
        if line1.startswith(">"):
            pass
        else:
            fh.write(line1+" " +line2+ "\n")
            count += 1

    fh.close()
    print("Saved %s adapters in adapters_combined.fa" % count)


def adapters_to_clean_ngs(filename):
    fh1 = open(filename, 'r')
    data1 = fh1.readlines()
    fh1.close()

    count = 0
    fh = open("adapters_ngs.txt", "w")
    for line in data1:
        line = line.strip().strip("\n")
        if line.startswith('>'):
            pass
        else:
            data = "adapter_%s\t%s\t0.5\t31\t10\t0\t0\n"% (count+1, line)
            fh.write(data)
            count+=1
    fh.close()



#adapters_to_clean_ngs("adapters_48_PCR-free_FWD.fa")
#if __name__ == "__main__":
#    import sys
#    args = sys.argv
#    adapters_files_to_list(args[1], args[2])


def adapter_removal_parser(filename):
    """Parses output of AdapterRemoval"""
    results = {}

    with open(filename, "r") as fin:
        lines = fin.readlines()
        for line in lines:
            if line.startswith("  --adapter"):
                lhs, rhs = line.split(":")
                name = lhs.strip().replace("-", "")
                sequence = rhs.strip()
                results[name] = sequence
    return results




class AdapterDB(object):
    """

    The name of the Fasta should be formatted as

        >Name|kraken:taxid|id

    where id is a number starting with 1000000
    This convention is adopted for now but may change. 

    """
    def __init__(self, filename=None):

        self.df = pd.DataFrame(columns=["name", "sequence", 
            "comment", "identifier", "filename"])

        if filename:
            self.load_fasta(filename)

    def load_all(self):
        from sequana.resources.data import adapters as dict_adapters
        from sequana import sequana_data
        for k,v in dict_adapters.items():
            self.load_fasta(sequana_data("data/%s" % v))

    def load_fasta(self, filename):
        from sequana.fasta import FastA
        adapters = FastA(filename)
        self.records = []
        for adapter in adapters:
            identifier = adapter.name.split("|")[2]
            record = {
                'name': adapter.name, 
                "sequence": adapter.sequence,
                "comment":adapter.comment,
                "identifier":identifier,
                "filename":filename}
            self.records.append(record)

        self.df = self.df.append(self.records)
        self.df.reset_index(drop=True, inplace=True)

        # check that identifiers should be unique
        if len(self.df) > len(self.df.identifier.unique()):
            print("Warn: there are duplicated identifiers in the adapters")

    def get_name(self, identifier):
        name =  self.df[self.df.identifier == str(identifier)].comment
        if len(name) == 1:
            name = list(name)[0]
        return name






