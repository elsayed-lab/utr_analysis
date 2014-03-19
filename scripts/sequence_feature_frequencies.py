#!/usr/bin/env python
"""
Keith Hughitt
2014/03/17

Finds reads containing either Poly(A) or Poly(T) tracts and attempts to guess
whether those sequences are a part of the genome, or come from a
poly-adenylation site.

The purpose of this script is to try and provide some basic guidelines for
where to look for poly-adenylation reads.
"""
import os
import glob
import re
from Bio import SeqIO,Seq

def main():
    # Select species and sample to query
    target = 'tcruzi' # 'lmajor'

    if target == 'lmajor':
        # Samples to query
        base_dir = "/cbcb/lab/nelsayed/raw_data/lminfectome"
        hpgl_id = "HPGL0075"  # procyclic (pathogen only)

        # output directory
        outdir = os.path.join('output', 'lmajor_hpgl0075')

        # Genome
        genome = os.path.join("/cbcb/lab/nelsayed/ref_data/lmajor_friedlin/genome",
                              "TriTrypDB-7.0_LmajorFriedlin_Genome.fasta")

        # L. major SL sequence and its reverse complement
        sl = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG"
        reverse_sl = "CAATAAAGTACAGAAACTGATACTTATATAGCGTTAGTT"
    else:
        # Samples to query
        base_dir = "/cbcb/lab/nelsayed/raw_data/tcruzir21"
        hpgl_id = "HPGL0250"  # trypomastigote (pathogen only)

        # output directory
        outdir = os.path.join('output', 'tcruzi_hpgl0250')

        # Genome
        genome = os.path.join("/cbcb/lab/nelsayed/ref_data/tcruzi_clbrener/genome/tc_esmer",
                              "TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta")

        # T. cruzi SL sequence and reverse complement
        sl = "AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG"
        reverse_sl = "CAATATAGTACAGAAACTGTATCAATAATAGCGTTAGTT"

    # RNA-Seq read filepaths
    reads = glob.glob(os.path.join(base_dir, hpgl_id, 'processed/*.fastq'))

    # regular expressions
    min_length = 16

    # spliced leader regular expressions
    sl_regex = '|'.join(["^" + sl[-x:] for x in range(min_length, len(sl) + 1)])
    reverse_sl_regex = '|'.join(
        [reverse_sl[-x:] + "$" for x in range(min_length, len(reverse_sl) + 1)]
    )

    search_patterns = {
        "polya_left": "^A{%d,}" % min_length,
        "polyt_left": "^T{%d,}" % min_length,
        "polya_right": "A{%d,}$" % min_length,
        "polyt_right": "T{%d,}$" % min_length,
        "sl_left": sl_regex,
        "rcsl_right": reverse_sl_regex
     }

    # count occurances of each sequence in reads
    for name,regex in search_patterns.items():
        print("Processing %s" % name)
        output_dir = os.path.join(outdir, name)
        count_seq_hits(re.compile(regex), reads, genome, output_dir)

def count_seq_hits(regex, reads, genome, outdir, max_reads=float('inf')):
    """Counts the number of occurances of a specified sequence in a collection
       of reads."""
    # load genome as a list of chromosome SeqRecords
    chromosomes = list(SeqIO.parse(genome, format='fasta'))

    # lists to keep track of different types of read matches
    with_feature = []
    with_feature_rc = []
    without_feature = []
    without_feature_rc = []
    no_match = []

    # Iterate through sample files
    for filepath in reads:
        for i, entry in enumerate(readfq(open(filepath)), 1):
            # stop once we have reached desired number of reads
            if i > max_reads:
                break

            # get read sequence and id
            read_id = entry[0]
            read = entry[1]

            # check for sequence pattern in read
            match = re.search(regex, read)

            # stop here if read doesn't contain sequence of interest
            if match is None:
                continue

            # otherwise check genome for read and trimmed read
            if (len(read) - match.end()) >= match.start():
                trimmed_read = read[match.end():]
            else:
                trimmed_read = read[:match.start()]

            # reverse complement
            rc = str(Seq.Seq(read).reverse_complement())
            trimmed_rc = str(Seq.Seq(trimmed_read).reverse_complement())

            # if found, see if it appears in the genome as-is, or when trimmed
            # to remove matched feature
            genome_match = False

            for chromosome in chromosomes:
                # full read
                if chromosome.seq.count(read) > 0:
                    with_feature.append(read_id)
                    genome_match = True
                # full read (reverse complement)
                elif chromosome.seq.count(rc) > 0:
                    with_feature_rc.append(read_id)
                    genome_match = True
                # trimmed read
                elif chromosome.seq.count(trimmed_read) > 0:
                    without_feature.append(read_id)
                    genome_match = True
                # trimmed read (reverse complement)
                elif chromosome.seq.count(trimmed_rc) > 0:
                    without_feature_rc.append(read_id)
                    genome_match = True

                # stop checking once match is found in genome
                if genome_match is True:
                    break

            # no match
            if not genome_match:
                no_match.append(read_id)

    # Save output
    if not os.path.exists(outdir):
        os.makedirs(outdir, mode=0o755)

    fp = open(os.path.join(outdir, 'full_read_matches.txt'), 'w')
    fp.write('\n'.join(with_feature) + '\n')
    fp.close()

    fp = open(os.path.join(outdir, 'full_read_reverse_matches.txt'), 'w')
    fp.write('\n'.join(with_feature_rc) + '\n')
    fp.close()

    fp = open(os.path.join(outdir, 'trimmed_read_matches.txt'), 'w')
    fp.write('\n'.join(without_feature) + '\n')
    fp.close()

    fp = open(os.path.join(outdir, 'trimmed_read_reverse_matches.txt'), 'w')
    fp.write('\n'.join(without_feature_rc) + '\n')
    fp.close()

    fp = open(os.path.join(outdir, 'no_matches.txt'), 'w')
    fp.write('\n'.join(no_match) + '\n')
    fp.close()

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
