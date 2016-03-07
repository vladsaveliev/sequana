import os

from sequoia import BAM


__all__ = ['bam_to_mapped_unmapped_fastq']

def bam_to_mapped_unmapped_fastq(filename, mode='pe'):
    """

    Given a BAM file, create FASTQ with R1/R2 reads mapped and unmapped.
    In the paired-end case, 4 files are created.


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

        for i,this in enumerate(bam):
            # check that this is paired indeed
            if this.is_paired is False:
                unpaired +=1
                # remove it ?
            else:
                # inpysam, seq is a string and qual a bytes....
                txt = b"@" + bytes(this.qname, "utf-8") + b"\t%s\n"
                txt += bytes(this.seq, "utf-8") + b"\n"
                txt += b"+\n"
                txt += this.qual + b"\n"

                if this.is_read1:
                    if this.is_unmapped:
                        R1_unmapped.write(txt % b"1:N:0:GTGAAA")
                        stats['R1_unmapped'] += 1
                    else:
                        R1_mapped.write(txt % b"1:N:0:GTGAAA")
                        stats['R1_mapped'] += 1
                elif this.is_read2:
                    if this.is_unmapped:
                        R2_unmapped.write(txt % b"2:N:0:GTGAAA")
                        stats['R2_unmapped'] += 1
                    else:
                        R2_mapped.write(txt % b"2:N:0:GTGAAA")
                        stats['R2_mapped'] += 1
                if this.is_duplicate:
                    print(this.qname)
        print(unpaired)
        print(stats)
        R1_mapped.close()
        R2_mapped.close()
        R1_unmapped.close()
        R2_unmapped.close()
    elif mode == "se":
        raise NotImplementedError
    else:
        raise ValueError("mode must be se or pe")

