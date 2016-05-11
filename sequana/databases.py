


"""
http://www.ebi.ac.uk/genomes/archaea.html

"""

import os
import glob

try:
    # py2
    from urllib.request import urlopen
except:
    from urllib import urlopen


from bioservices import EUtils, ENA
from easydev import AttrDict


class EUtilsTools(object):

    def __init__(self):
        self.eutils = EUtils()

    def accession_to_info(self, ids):
        print('Fetching using EUtils')
        res = self.eutils.EFetch(db="nuccore",id=ids, 
                rettype="docsum", retmode="dict")

        res = res['eSummaryResult']['DocSum']

        # if one id provided, it will be a dict, otherwise a list of dicts
        try:
            res[0]
        except:
            res = [res]

        # now we can loop over all identifiers
        records = {}
        accessions = [x.strip() for x in ids.split(',')]

        for i, entry in enumerate(res):
            # first, save the acc number
            accession = entry['Id']
            # then various info
            items = entry['Item']
            identifier = [x for x in items if x['@Name'] == "Extra"][0]['#text']
            if "||" in identifier:
                # strip content after || 
                identifier = identifier.split("||")[0] 

            title = [x for x in items if x['@Name'] == "Title"][0]['#text']
            taxid = [x for x in items if x['@Name'] == "Title"][0]['#text']
            gi = [x for x in items if x['@Name'] == "Gi"][0]['#text']
            record = {
                "taxid": taxid,
                'accession': accession, 
                "identifier": identifier,
                'gi': gi,
                'comment': title
                }

            records[accessions[i]] = AttrDict(**record)
        return records




class ENADownload(object):
    """

    Downloader facility to retrieve genome fasta files from ENA


    """
    def __init__(self):
        self.convert_enaacc_to_gi = True
        self.eutils = EUtilsTools()

    def download_fasta(self, filelist, output_dir=None):
        """

        :param filelist: a name to find on the ENA web server OR the 
            name of an accession number.

        """
        from bioservices import ENA
        if filelist.endswith(".txt"):
            print("Downloading list from http://www.ebi.ac.uk/genomes/%s" % filelist)
            data = urlopen("http://www.ebi.ac.uk/genomes/%s" % filelist).readlines()
            identifiers = [x.strip().decode() for x in data]
        else:
            identifiers = [filelist]

        ena = ENA()

        if output_dir is None:
            output_dir = "."
        else:
            try: os.mkdir(output_dir)
            except:pass

        from easydev import Progress
        N = len(identifiers)
        pb = Progress(N)
        for i, identifier in enumerate(identifiers):
            filenames = glob.glob(output_dir + os.sep + "ENA_%s*" % identifier)

            if len(filenames) >= 1:
                # no need to fetch and save the data it looks like...
                continue

            # download data from ENA
            data = ena.get_data(identifier, "fasta")

            # Split header and Fasta
            header, others = data.decode().split("\n", 1)

            name  = header.strip(">").split(" ")[0]
            db, id_, acc = name.split("|")

            # Here, we will fetch info from NCBI every time. This is a pity but
            # no choice.
            header = self.switch_header_to_gi(acc)

            # Save to local file
            filename = "%s_%s.fasta" % (db, acc.split(".")[0])
            if output_dir:
                filename = output_dir + os.sep + filename

            with open(filename, "w") as fout:
                fout.write(header+"\n"+others)
            pb.animate(i+1)

    def switch_header_to_gi(self, acc):
        res = self.eutils.accession_to_info(acc)[acc]
        return ">"+res['identifier']+" " + res['comment']

    def download_viruses(self):
        """4700 organisms (may 2016)"""
        self.download_fasta("virus.txt", "Viruses")

    def download_plasmids(self):
        """2800 organisms (may 2016) """
        self.download_fasta("plasmid.txt", "Plasmids")

    def download_phage(self):
        """2400 organisms (may 2016) """
        self.download_fasta("phage.txt", "Phages")

    def download_archaealvirus(self):
        """2400 organisms (may 2016) """
        self.download_fasta("archaealvirus.txt", "ArchaeaViruses")

    def download_archaea(self):
        """2400 organisms (may 2016) """
        self.download_fasta("archaea.txt", "Archaea")

    def download_bacteria(self):
        """2400 organisms (may 2016) 

        THis is the longest and takes about 20mins
        """
        self.download_fasta("bacteria.txt", "Bacteria")

    def download_escheveria_coli(self):
        # Included in bacteria
        pass

    def download_accession(self, acc):
        self.download_fasta(acc, "Custom")




