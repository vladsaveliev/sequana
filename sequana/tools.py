import os
import string
import glob

from sequana import BAM


__all__ = ['bam_to_mapped_unmapped_fastq']


try:
    _translate = string.maketrans('ACGTacgt', 'TGCAtgca')
    def bytes(x,dummy):
        return x
except:
    _translate = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')


def reverse_complement(seq):
    # use hastag but for now, do it manually
    return seq.translate(_translate)[::-1]



def bam_to_mapped_unmapped_fastq(filename, output_directory=None):
    """Create mapped and unmapped fastq files from a BAM file

    Given a BAM file, create FASTQ with R1/R2 reads mapped and unmapped.
    In the paired-end case, 4 files are created.

    The interest of this function is that it does not create intermediate files
    limiting IO in the process.

    Reads that are not paired or not "proper pair" (neither flag 4 nor flag 8)
    are ignored

    Secondary alignment (flag 256) are dropped so as to remove any ambiguous alignments and
    keep the number of final reads equal to the initial number of reads

    If R1 is mapped or R2 is mapped then the reads are considered mapped. If
    both R1 and R2 are unmapped, then reads are unmapped.

    Note about chimeric alginment: one is the representative and the other is
    the supplementary. This flag is not used in this function. Note also that
    chimeric alignment have same QNAME and flag 4 and 8

    .. todo:: comments is currently harcoded and should be removed. comments
        from the fastq qname are not stored in BAM.
    """
    bam = BAM(filename)
    # figure out if this is paired or unpaired

    newname, ext = os.path.splitext(filename)

    import collections
    stats = collections.defaultdict(int)
    stats['R1_unmapped'] = 0
    stats['R1_mapped'] = 0


    # figure out where to save the file 
    if output_directory is None:
        pass
    else:
        assert isinstance(filename, str)
        from sequana.snaketools import FileFactory
        ff = FileFactory(filename)
        newname = output_directory + os.sep + ff.filenames[0]

    R1_mapped = open(newname + "_R1.mapped.fastq", "wb")
    R1_unmapped = open(newname + "_R1.unmapped.fastq", "wb")
    stats['duplicated'] = 0
    stats['unpaired'] = 0

    unpaired = 0

    # if paired, let open other files
    if bam.is_paired:
        stats['mode'] = "pe"
        stats['R2_unmapped'] = 0
        stats['R2_mapped'] = 0
        R2_mapped = open(newname + "_R2.mapped.fastq", "wb")
        R2_unmapped = open(newname + "_R2.unmapped.fastq", "wb")
    else:
        stats['mode'] = "se"


    # loop through the BAM (make sure it is rewind)
    bam.reset()

    from easydev import Progress
    pb = Progress(len(bam))
    for i, this in enumerate(bam):
        if this.flag & 256:
            # Unmapped reads are in the BAM file but have no valid assigned
            # position (N.B., they may have an assigned position, but it should be ignored).
            # It's typically the case that a number of reads can't be aligned, due to things
            # like sequencing errors, imperfect matches between the DNA sequenced and the
            # reference, random e. coli or other contamination, etc..
            # A secondary alignment occurs when a given read could align reasonably well to
            # more than one place. One of the possible reported alignments is termed "primary"
            # and the others will be marked as "secondary".
            stats['secondary'] +=1 
            if this.is_paired is False:
                stats['unpaired'] += 1
        else:
            # inpysam, seq is a string and qual a bytes....
            if this.is_reverse is True:
                #print("reversed complement entry %s" % i)
                txt = b"@" + bytes(this.qname, "utf-8") + b"\n"
                revcomp = reverse_complement(this.seq)
                txt += bytes(revcomp, "utf-8") + b"\n"
                txt += b"+\n"
                txt += bytes(this.qual[::-1], 'utf-8') + b"\n"
            else:
                txt = b"@" + bytes(this.qname, "utf-8") + b"\n"
                txt += bytes(this.seq, "utf-8") + b"\n"
                txt += b"+\n"
                txt += bytes(this.qual,"utf-8") + b"\n"

            # Here, we must be careful as to keep the pairs. So if R1 is mapped
            # but R2 is unmapped (or the inverse), then the pair is mapped
            if this.is_read1:
                if this.is_unmapped and this.mate_is_unmapped:
                    R1_unmapped.write(txt)
                    stats['R1_unmapped'] += 1
                else:
                    R1_mapped.write(txt)
                    stats['R1_mapped'] += 1
            elif this.is_read2:
                if this.is_unmapped and this.mate_is_unmapped:
                    R2_unmapped.write(txt)
                    stats['R2_unmapped'] += 1
                else:
                    R2_mapped.write(txt)
                    stats['R2_mapped'] += 1
            else:
                # This should be a single read
                #assert self.is_paired is False
                stats['unpaired'] += 1
                if this.is_unmapped:
                    R1_unmapped.write(txt)
                    stats['R1_unmapped'] += 1
                else:
                    R1_mapped.write(txt)
                    stats['R1_mapped'] += 1

            if this.is_duplicate:
                stats['duplicated'] += 1
        pb.animate(i+1)

    if bam.is_paired:
        R2_mapped.close()
        R2_unmapped.close()
    print("\nNumber of entries in the BAM: %s" % str(i+1))
    R1_mapped.close()
    R1_unmapped.close()
    return stats
bam2fastq = bam_to_mapped_unmapped_fastq




def bam_get_paired_distance(filename):
    """


    return position start and end of the paired-end reads that were mapped
    (both)

    """

    b = BAM(filename)
    distances = []

    for fragment in b:
        if fragment.is_unmapped is False and fragment.mate_is_unmapped is False \
            and fragment.is_read1:

            # get the mate:
            mate = next(b)


            if fragment.is_reverse:
                position2 = fragment.reference_end
                position1 = mate.reference_start
                mode = 1
            elif mate.is_reverse:
                position1 = fragment.reference_start
                position2 = mate.reference_end
                mode = 2
            else: # if both are not reversed, what does that mean.
                # On Hm2, this is the case for 4 pairs out of 1622
                # This seems to be a special case for fragment ends exactly 
                # at the end of the reference and mate starts exactly at 
                # the beginnin with a length less than 100
                print(fragment.reference_start, fragment.reference_end)
                print(mate.reference_start, mate.reference_end)

            distances.append((position1, position2, mode))

    return distances


class FastqFactory(dict):
    """List a set of files and identify projects (paired or not)

    Fastq files may be provided within a single directory and one
    would need to identify the different paired-files (or not).

    This factory is a dictionary and project names with be the keys and the
    values contain the list of files (either one or two).


    """
    def __init__(self, pattern, compression_ext="gz"):
        """.. rubric:: Constructor

        :param str pattern: a valid pattern (e.g., '*fastq*') to a list of 
            FastQ files
        :param str compression_ext: the possible compression extension. 
        """
        self._compext = compression_ext
       
        if isinstance(pattern, list):
            self.filenames = pattern
        elif isinstance(pattern, str):
            self.filenames = glob.glob(pattern)


        if len(self.filenames) ==0:
            raise ValueError("No files found")
        # checks that they ends in either .gz or .fastq
        for filename in self.filenames:
            if filename.endswith('.'+self._compext) or filename.endswith(".fastq"):
                pass
            else:
                raise ValueError("One of the filename does not end "
                "in .fastq or .fastq.%s (%s)" % (self._compext, filename))

        from .snaketools import FileFactory
        self.ff = FileFactory(pattern)

        for project in self.project_names:
            self[project] = [os.sep.join([path]) for base, path in 
                zip(self.ff.basenames, self.ff._glob) if base.startswith(project)]


    def _get_project_names(self):
        # strip the .fastq.gz or .fastq if any and actually even
        # .blabla.fastq.gz
        project_names = self.ff.filenames

        # all projects should contains either one fastq (R1) or a pair (R1/R2)
        # The project name is the string before _R1
        for name in project_names:
            if "_R1" not in name and "_R2" not in name:
                raise ValueError("Only Fastq with the tag _R1_ or _R2_ are handled")

        project_names = [x.split("_R1")[0] for x in project_names]
        project_names = [x.split("_R2")[0] for x in project_names]

        project_names = sorted(list(set(project_names)))
        return project_names
    project_names = property(_get_project_names) 

















