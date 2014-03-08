#!/usr/bin/env python
"""
Keith Hughitt
2014/01/29

Counts the number of occurrences of each 10+ nt suffix of the T. cruzi spliced
leader sequence in the genome. This should provide us with some sense of how
likely a non-left-aligned SL subsequence is to be part of the SL (with some
artifact sequence to the left), or a similar but unrelated sequence.
"""
import os
from Bio import SeqIO

sl="AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG"

# load genome
filepath = os.path.join("/cbcb/lab/nelsayed/ref_data/tcruzi_clbrener/genome/tc_esmer",
                        "TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta")
genome = list(SeqIO.parse(filepath, format='fasta'))

counts = []

# iterate through each 10+ nt suffix of the SL and count the number of
# occurrences in the genome
suffixes = [sl[-n:] for n in range(8, len(sl) + 1)]

for suffix in suffixes:
    num_occurrences = 0
    for chromosome in genome:
        num_occurrences += chromosome.seq.count(suffix)

    counts.append(num_occurrences)

# print results
print("Length | Suffix: Number of Occurrences")
for i, suffix in enumerate(suffixes):
    print("%2d | %s:%d" % (len(suffix), suffix.rjust(40), counts[i]))

