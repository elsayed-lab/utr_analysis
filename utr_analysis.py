#!/bin/env python
# -*- coding: utf-8 -*-
"""
Utranslated Region (UTR) analysis
Keith Hughitt
2013/12/15

The purpose of this script is the scan a collection of RNA-Seq reads and
determine the location of the 5'UTR splice acceptor sites or 3'UTR
poly-adenylation sites.

A reference genome and CDS coordinates are required as input and will be used
to map reads back after removing the spliced leader in order to determine the
SL acceptor site.

NOTE:

Initially, development of this script will focus on SL site determination. Once
this has been properly implemented, support for poly-A analysis will be added.

Testing:

utr_analysis.py \
    -i $RAW/tcruzir21/HPGL0258/processed/*.filtered.fastq \ 
    -f $REF/tcruzi_clbrener/genome/tc_esmer/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta \
    -g $REF/tcruzi_clbrener/annotation/tc_esmer/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like.gff      \
    -s AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG output.csv

"""
import os
import sys
import re
import glob
import textwrap
import argparse

def main():
    """Main"""
    args = parse_input()
    matches = parse_reads(args)

    # save output to a file for now
    with open(args.output, 'w') as fp:
        fp.writelines(matches)

def parse_reads(args):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing the feature of interest (SL or polyA tail).
    """
    matches = []

    # limit to matches of size min_length or greater
    s = spliced_leader[-args.min_length:]

    # compile regex allowing a single mismatch;
    # we will filter based on number of mistmatches again later on a larger
    # portion of the read using the value specified at run-time
    regex = re.compile('|'.join('^%s.%s' % (s[:i], s[i+1:]) for i in range(len(s))))

    for filepath in glob.glob(args.input_reads):
        # open fastq file
        fastq = open(filepath)
        print("Processing %s" % filepath)

        # filter by length
        for i, read in enumerate(readfq(fastq)):
            seq = read[1][:len(args.spliced_leader)]

            # check for match
            if re.search(regex, seq) is not None:
                matches.append([i, seq])

    # filtered sequeneces
    return matches


def parse_input():
    """
    Parses script input and returns values.
    """

    # Usage example
    usage_examples=textwrap.dedent("""\
    Usage Example:
    --------------
    ./utr_analysis.py                                              \\
        -i $RAW/tcruzir21/*/processed/*.filtered.fastq             \\
        -s AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG                 \\
        -f TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta \\
        -g TrypDB-6.0_TcruziCLBrenerEsmeraldo-like.gff             \\
        -m 12 -n 3 output.csv
    """)

    # Create ArgumentParser instance
    parser = argparse.ArgumentParser(
        description='Spliced Leader and poly-adenylation site analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=usage_examples
    )

    # Add arguments
    parser.add_argument('-i', '--input-reads', required=True,
                        help='RNA-Seq FASTQ glob string')
    parser.add_argument('-f', '--fasta-genome', dest='fasta', required=True,
                        help='Genome sequence FASTA filepath')
    parser.add_argument('-g', '--gff-annotation', dest='gff', required=True,
                        help='Genome annotation GFF')
    parser.add_argument('-s', '--sl-sequence', dest='spliced_leader', 
                        help='Spliced leader DNA sequence')
    parser.add_argument('-m', '--min-length', default=10, type=int,
                        help='Minimum length of SL match (default=10)')
    parser.add_argument('-n', '--num-mismatches', default=2, type=int,
                        help='Number of mismatches to allow (default=2)')
    parser.add_argument('output', metavar='OUTPUT FILENAME',
                        help='Filepath to write CSV output to.')

    # Parse arguments
    args = parser.parse_args()

    # Replace any environmental variables and return args
    args.fasta = os.path.expandvars(args.fasta)
    args.input_reads = os.path.expandvars(args.input_reads)
    args.gff = os.path.expandvars(args.gff)
    args.output = os.path.expandvars(args.output)

    # @TODO Validate input

    return args

def readfq(fp):
    """
    Loads a fastq file into memory
    https://github.com/lh3/readfq
    """
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
    sys.exit(main())
