


class CutAdaptParser(object):

    def __init__(self, filename):
        self.filename = filename
        with open(self.filename) as fin:
            data = fin.read()
            if "Total read pairs processed" in data:
                self.mode = "pe"
            else:
                self.mode = 'se'
        print("Mode: %s" % self.mode)

    def _get_data_tobefound(self):
        tobefound = []
        if self.mode == 'se':
            tobefound.append(('total_reads', 'Total reads processed:'))
            tobefound.append(('reads_with_adapters', 'Reads with adapters:'))
            tobefound.append(('reads_too_short', 'Reads that were too short:'))
            tobefound.append(('reads_too_long', 'Reads that were too long:'))
            tobefound.append(('reads_kept', 'Reads written (passing filters):'))
            #tobefound.append(('total_basepairs_processed', 'Total basepairs processed:'))
            #tobefound.append(('total_basepairs_filtred', 'Total written (filtered):'))
        else:
            tobefound.append(('total_reads', 'Total read pairs processed:'))
            tobefound.append(('reads_too_short', 'Pairs that were too short:'))
            tobefound.append(('reads_kept', 'Pairs written (passing filters):'))

        return tobefound

    def parse(self):
        d = {}

        tobefound = self._get_data_tobefound()

        adapters = []

        with open(self.filename, 'r') as fin:
            # not tool large so let us read everything
            data = fin.readlines()

            # some metadata to extract
            for this in tobefound:
                key, pattern = this
                found = [line for line in data if line.startswith(pattern)]
                if len(found) == 0:
                    print("%s not found)" % pattern)
                elif len(found) == 1:
                    d[key] = found[0].split(":", 1)[1].strip()

            adapters = []
            dd = {}
            for this in data:
                if this.startswith("=== Adapter"):
                    name = this.split("=== Adapter")[1].split(" ===")[0].strip()
                    dd['name'] = name
                    continue
                if this.startswith('Sequence:'):
                    info = this.split("Sequence:", 1)[1].strip()
                    info = info.split(";")
                    dd['info'] = {
                        'Sequence': info[0].strip(),
                        'Type': info[1].split(':',1)[1].strip(),
                        'Length': info[2].split(':',1)[1].strip(),
                        'Trimmed': info[3].split(':',1)[1].strip()
                    }
                    
                    adapters.append(dd.copy())

        return d, adapters

    def table(self):
        pass
