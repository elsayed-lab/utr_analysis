#!/usr/bin/env python
"""
Keith Hughitt
2014/02/17

Counts the number of occurrences of each 10+ nt suffix of the L. major spliced
leader (SL) sequence, SL sequence, and Poly(A) and Poly(T) tracts within a
randomly chosen L. major RNA-Seq sample.

The purpose of this script is provide a sense for the possible distribution
of actual UTR feature sequences of interest, where they can be most easily
observed, and how many hits for the sequence of interest (or possible unrelated
false hits) occur at each subsequence length.
"""
import os
import glob
import re

def main():
    # Samples to query
    #base_dir = "/cbcb/lab/nelsayed/raw_data/lminfectome"
    #hpgl_id = "HPGL0075"  # procyclic (pathogen only)

    # L. major SL sequence and its reverse complement
    #sl = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG"
    #reverse_sl = "CAATAAAGTACAGAAACTGATACTTATATAGCGTTAGTT"

    # Samples to query
    base_dir = "/cbcb/lab/nelsayed/raw_data/tcruzir21"
    hpgl_id = "HPGL0250"  # trypomastigote (pathogen only)

    # T. cruzi SL sequence and reverse complement
    sl = "AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG"
    reverse_sl = "CAATATAGTACAGAAACTGTATCAATAATAGCGTTAGTT"

    # Min suffix length
    min_length = 10

    # Number of reads to query (can set a limit to speed things up)
    max_reads = float('inf')

    # list to keep track of counts for each subsequence length
    starting_counter = {i:0 for i in range(min_length, len(sl) + 1)}

    counts = {
        'sl': {
            'left': starting_counter.copy(),
            'right': starting_counter.copy(),
            'any': starting_counter.copy()
        },
        'reverse_sl': {
            'left': starting_counter.copy(),
            'right': starting_counter.copy(),
            'any': starting_counter.copy()
        },
        'polya': {
            'left': starting_counter.copy(),
            'right': starting_counter.copy(),
            'any': starting_counter.copy()
        },
        'polyt': {
            'left': starting_counter.copy(),
            'right': starting_counter.copy(),
            'any': starting_counter.copy()
        }
    }

    # R1 and R2 filepaths
    files = glob.glob(os.path.join(base_dir, hpgl_id, 'processed/*.fastq'))

    # iterate through each 10+ nt suffix of the SL and count the number of
    # occurrences in the genome
    sl_suffixes = [sl[-n:] for n in range(min_length, len(sl) + 1)]
    reverse_prefixes = [reverse_sl[:n] for n in range(min_length, len(sl) + 1)]

    for i, sl_suffix in enumerate(sl_suffixes, min_length):
        # Spliced leader regexes
        features = {
            'sl': sl_suffix,
            'reverse_sl': reverse_prefixes[i - min_length],
            'polya': "A" * i,
            'polyt': "T" * i
        }

        # Iterate through sample files
        for filepath in files:
            for j, read in enumerate(readfq(open(filepath)), 1):
                # stop once we have reached desired number of reads
                if j > max_reads:
                    break

                # Check for features at each possible position in read
                for feature in features:
                    # first see if it exists at all in read
                    if features[feature] in read[1]:
                        counts[feature]['any'][i] += 1

                        # if so, check at ends as well
                        if read[1].startswith(features[feature]):
                            counts[feature]['left'][i] += 1
                        if read[1].endswith(features[feature]):
                            counts[feature]['right'][i] += 1

    # Print header
    header_fields = ['length']
    for feature in counts:
        for pos in counts[feature]:
            header_fields.append("_".join([feature, pos]))
    print(",".join(header_fields))

    # Print rows
    for i in range(min_length, len(sl) + 1):
        row = [str(i)]
        for feature in counts:
            for pos in counts[feature]:
                row.append(str(counts[feature][pos][i]))
        print(",".join(row))

# FASTQ parser
# source: https://github.com/lh3/readfq
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

if __name__ == "__main__":
    main()
