import os
import string

from sequoia import BAM


__all__ = ['bam_to_mapped_unmapped_fastq']


try:
    _translate = string.maketrans('ACGTacgt', 'TGCAtgca')
except:
    _translate = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')


def reverse_complement(seq):
    # use hastag but for now, do it manually
    return seq.translate(_translate)[::-1]

def bam_to_mapped_unmapped_fastq(filename, mode='pe'):
    """

    Given a BAM file, create FASTQ with R1/R2 reads mapped and unmapped.
    In the paired-end case, 4 files are created.

    Reads that are not paired or not "proper pair" (neither flag 4 nor flag 8)
    are ignored


    This is equivalent to the following script (for example)::

        # map alignments (4 is flag for unmapped segment and -F is do not
        # include them)
        samtools view -F 4 filename.bam -o mapped.sam
        # unmapped segments
        samtools view -F 4 filename.bam -o unmapped.sam
        awk '{print $1}' unmapped.sam | uniq  > unmapped.bed
        awk '{print $1}' mapped.sam | uniq > mapped.bed

        # create the fastq files based on the fastq and bed files
        seqtk subseq {input.R1} unmapped.bed > R1.unmapped.fastq
        seqtk subseq {input.R1} mapped.bed > R1.mapped.fastq
        seqtk subseq {input.R2} unmapped.bed > R2.unmapped.fastq
        seqtk subseq {input.R2} mapped.bed > R2.mapped.fastq

    This function does not create intermediate files.


    """
    bam = BAM(filename)

    newname, ext = os.path.splitext(filename)

    if mode == "pe":
        import collections
        stats = collections.defaultdict(int)
        unpaired = 0

        R1_mapped = open(newname + "_R1.mapped.fastq", "bw")
        R2_mapped = open(newname + "_R2.mapped.fastq", "bw")
        R1_unmapped = open(newname + "_R1.unmapped.fastq", "bw")
        R2_unmapped = open(newname + "_R2.unmapped.fastq", "bw")

        

        from easydev import Progress
        pb = Progress(len(bam))
        for i, this in enumerate(bam):
            # check that this is paired indeed
            if this.is_paired is False:
                stats['unpaired'] +=1
                # What to do in such case ?
            elif this.flag & 256:
                # Unmapped reads are in the BAM file but have no valid assigned
                # position (N.B., they may have an assigned position, but it should be ignored).
                # It's typically the case that a number of reads can't be aligned, due to things
                # like sequencing errors, imperfect matches between the DNA sequenced and the
                # reference, random e. coli or other contamination, etc..
                # A secondary alignment occurs when a given read could align reasonably well to
                # more than one place. One of the possible reported alignments is termed "primary"
                # and the others will be marked as "secondary".

                stats['secondary'] +=1 
            else:
                # inpysam, seq is a string and qual a bytes....
                if this.is_reverse is True:
                    #print("reversed complement entry %s" % i)
                    txt = b"@" + bytes(this.qname, "utf-8") + b"\t%s\n"
                    revcomp = reverse_complement(this.seq)
                    txt += bytes(revcomp, "utf-8") + b"\n"
                    txt += b"+\n"
                    txt += this.qual[::-1] + b"\n"
                else:
                    txt = b"@" + bytes(this.qname, "utf-8") + b"\t%s\n"
                    txt += bytes(this.seq, "utf-8") + b"\n"
                    txt += b"+\n"
                    txt += this.qual + b"\n"

                # Here, we must take keep the pairs so if R1 is mapped
                # but R2 is unmapped (or the inverse), then the pair is
                # unmapped
                if this.is_read1:
                    if this.is_unmapped or this.mate_is_unmapped:
                        R1_unmapped.write(txt % b"1:N:0:GTGAAA")
                        stats['R1_unmapped'] += 1
                    else:
                        R1_mapped.write(txt % b"1:N:0:GTGAAA")
                        stats['R1_mapped'] += 1
                elif this.is_read2:
                    if this.is_unmapped or this.mate_is_unmapped:
                        R2_unmapped.write(txt % b"2:N:0:GTGAAA")
                        stats['R2_unmapped'] += 1
                    else:
                        R2_mapped.write(txt % b"2:N:0:GTGAAA")
                        stats['R2_mapped'] += 1
                else:
                    #pass
                    print("Neither R1 or R2")
                if this.is_duplicate:
                    stats['duplicated'] += 1
            pb.animate(i+1)
        print('')
        print("Number of entries in the BAM: %s" % str(i+1))
        R1_mapped.close()
        R2_mapped.close()
        R1_unmapped.close()
        R2_unmapped.close()
        return stats
    elif mode == "se":
        raise NotImplementedError
    else:
        raise ValueError("mode must be se or pe")


