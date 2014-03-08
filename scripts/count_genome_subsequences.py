#!/usr/bin/env python
"""
Keith Hughitt
2014/03/07

Counts the number of occurrences of each 10+ nt suffix of the T. cruzi spliced
leader sequence in the genome, as well as the number of occurences of Poly(A),
Poly(T), Poly(G), and random sequences of similar lengths.

This should provide us with some sense of how likely reads with similar
sequence fragments are to have come from the genome rather than a spliced
leader or poly-A addition.
"""
import os
import random
from Bio import SeqIO

# load genome
filepath = os.path.join("/cbcb/lab/nelsayed/ref_data/tcruzi_clbrener/genome/tc_esmer",
                        "TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta")
genome = list(SeqIO.parse(filepath, format='fasta'))

# Spliced leader sequence (T. cruzi)
sl="AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG"
counts_sl = []

# iterate through each 10+ nt suffix of the SL and count the number of
# occurrences in the genome
suffixes = [sl[-n:] for n in range(10, len(sl) + 1)]

for suffix in suffixes:
    num_occurrences = 0
    for chromosome in genome:
        num_occurrences += chromosome.seq.count(suffix)

    counts_sl.append(num_occurrences)

counts_sl = counts_sl + [0] * 30 # always zero after length 39...

# Poly(A) occurences
polya = ["A" * i for i in range(10, 56)]
counts_a = []

for sequence in polya:
    num_occurrences = 0
    for chromosome in genome:
        num_occurrences += chromosome.seq.count(sequence)

    counts_a.append(num_occurrences)

# Poly(T) occurences
counts_t = []
polyt = ["T" * i for i in range(10, 56)]

for sequence in polyt:
    num_occurrences = 0
    for chromosome in genome:
        num_occurrences += chromosome.seq.count(sequence)

    counts_t.append(num_occurrences)

# Poly(G) occurences
counts_g = []
polyg = ["G" * i for i in range(10, 56)]

for sequence in polyg:
    num_occurrences = 0
    for chromosome in genome:
        num_occurrences += chromosome.seq.count(sequence)

    counts_g.append(num_occurrences)

# for comparison, look at occurences of poly-G tracts and random substrings
# as well
alphabet = ["A", "G", "C", "T"]
random_seq = "".join([alphabet[random.randint(0,3)] for x in range(56)])
print("Random sequence: %s" % random_seq)

counts_rand = []
random_subseqs = [random_seq[:i] for i in range(10, 56)]

for sequence in random_subseqs:
    num_occurrences = 0
    for chromosome in genome:
        num_occurrences += chromosome.seq.count(sequence)

    counts_rand.append(num_occurrences)

# print results
print("seq_length,matches_sl,matches_a,matches_t,matches_g,matches_rand")
for i, count in enumerate(counts_t):
    print("%2d,%d,%d,%d,%d,%d" % (len(polyt[i]), counts_sl[i], 
                                  counts_a[i], counts_t[i],
                                  counts_g[i], counts_rand[i]))

