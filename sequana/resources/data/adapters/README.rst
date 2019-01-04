The adapters are stored as FASTA files compatible with cutadapt.
The files are exemplars of what can be used and are not universal. 

They are provided in 3 different versions: standard, reverse and reverse
complement; these 3 files contain the same information.



Summary of what's included
--------------------------

- Nextera i5 Index 2 including all S/N indexes in adapters_Nextera
- Nextera i7 Index 1 including only N indexes in adapters_Nextera
- TruSeq Single Indexes in adapters_TruSeq
- TruSeq Small RNA in adapters_Small



adapters_NEBNext2_fwd.fa
adapters_Rubicon_fwd.fa
adapters_SMARTer_fwd.fa
adapters_NEBNext_fwd.fa
adapters_PCRFree_fwd.fa


Nextera DNA Indexes
--------------------
contains only Index 2 (i5) for now.

Transposase
~~~~~~~~~~~~~~~

The transposase adapters are used for Nextera tagmentation and are removed

Read 1::

    5′ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG

Read 2::

    5′ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

Indexes
~~~~~~~~~
The Nextera file contains i5 Bases for Sample Sheet MiSeq, HiSeq 2000/2500
as quoted in page 13 of http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf

known as S/E/N

N50X Nextera DNA, S50X Nextera XT and E50X Nextera Enrichment and Nextera Rapid Capture.

The series S/E/N have actually the same index.


PCRFree
----------
Bar code are 6-bases long
http://www.biooscientific.com/Portals/0/IEM/Bioo-Scientific-PCR-Free-Barcode-Indices-v1-1-15.pdf



16S
-------
http://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf

The dual indexing strategy uses two 8 base indices, Index1 (i7) adjacent to the P7 sequence, and Index2 (i5) adjacent to the P5 sequence. Dual indexing is enabled by adding a unique Index 1 (i7) and Index 2 (i5) to each  sample.

The  96 sample Nextera XT Index Kit (FC‐131- 1002) use12 different Index1 (i7) adapters (N701-N712) and 8 different Index2 (i5) adapters (S501-S508).

The 24 sample Nextera XT Index Kit (FC‐131-1001) uses 6 different Index 1 (i7) adapters (N701-N706) and 4 different Index 2 (i5) adapters (S501-S504). 

In the Index adapter name, the N or S refers
to Nextera XT sample preparation, 7 or 5 refers to Index1 (i7) or Index2 (i5), respectively.

The 01-12 refers to the Index number.

============ ========= ============ =======================
Index 1 (i7) Sequence, Index 2 (i5) Sequence
============ ========= ============ =======================
N701         TAAGGCGA  S501         TAGATCGC
N702         CGTACTAG  S502         CTCTCTAT
N703         AGGCAGAA  S503         TATCCTCT
N704         TCCTGAGC  S504         AGAGTAGA
N705         GGACTCCT  S505         GTAAGGAG
N706         TAGGCATG  S506         ACTGCATA
N707         CTCTCTAC  S507         AAGGAGTA
N708         CAGAGAGG  S508         CTAAGCCT
N709         GCTACGCT
N710         CGAGGCTG
N711         AAGAGGCA
N712         GTAGAGGA
============ ========= ============ =======================


TruSeq CD Indexes 
------------------------------------
:name: TruSeqCD

Combinatorial dual (CD) index adapters (formerly TruSeq HT). Extracted from IEM
using sequana.iem.IEM class.



TruSeq Single Indexes
----------------------
**current name is TruSeq**. Not to be confused with double indexing.


from TruSeq **LT** Kits and TruSeq v1/v2 Kits for DNA and RNA.

The following sequences are used for adapter trimming.

Read 1::

    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

Read 2::

    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

TruSeq UniversalAdapter::

    5′ AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT

Index adapter sequences are six bases. The index numbering is not sequential, 
so indexes 17, 24, and 26 are skipped. Additionally, the bases
preceding each index adapter sequence are the same, but the two bases following
the index adapter sequence can vary




TruSeq Small RNA
-----------------
from Illumina TruSeq small RNA Kits.

Based on document 10000002694 v09

!! in the document above, six-bases long index adapters are reverse
complemented. !!

Reference are eg RPI5 for index 5.


SMARTer
-------
from SMARTer Stranded RNA-Seq Kit


NEBNext
------------

Contains the single and dual-indexing


Universal found in manualE7335_index_primers_set1 of illumina.














