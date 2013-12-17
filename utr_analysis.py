#!/bin/env python2
# -*- coding: utf-8 -*-
"""
Utranslated Region (UTR) analysis
Keith Hughitt
2013/12/15

Overview
--------
The purpose of this script is the scan a collection of RNA-Seq reads and
determine the location of the 5'UTR splice acceptor sites or 3'UTR
poly-adenylation sites.

A reference genome and CDS coordinates are required as input and will be used
to map reads back after removing the spliced leader in order to determine the
SL acceptor site.

NOTE
----
Initially, development of this script will focus on SL site determination. Once
this has been properly implemented, support for poly-A analysis will be added.

Testing
-------

utr_analysis.py \
    -i "$RAW/tcruzir21/HPGL0258/processed/*.filtered.fastq" \ 
    -f $REF/tcruzi_clbrener/genome/tc_esmer/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta \
    -g $REF/tcruzi_clbrener/annotation/tc_esmer/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like.gff      \
    -s AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG output.csv
"""
import os
import sys
import re
import glob
import datetime
import textwrap
import argparse
from ruffus import *

def parse_input():
    """
    Parses script input and returns values.
    """
    # Usage example
    usage_examples=textwrap.dedent("""\
    Usage Example:
    --------------
    ./utr_analysis.py                                              \\
        -i "$RAW/tcruzir21/*/processed/*.filtered.fastq"           \\
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
    parser.add_argument('-a', '--author', help='Author contact name', 
                        default='')
    parser.add_argument('-e', '--email', help='Author contact email address',
                        default='')
    parser.add_argument('output', metavar='OUTPUT FILENAME',
                        help='Filepath to write CSV output to.')

    # Parse arguments
    args = parser.parse_args()

    # Replace any environmental variables and return args
    args.fasta = os.path.expandvars(args.fasta)
    args.input_reads = os.path.expandvars(args.input_reads)
    args.gff = os.path.expandvars(args.gff)
    args.output = os.path.expandvars(args.output)

    # set defaults for author if none is specified
    if args.author is "":
        args.author = os.getlogin()

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

def create_header_comment(filename, description, author, email):
    """
    Creates a header comment to be appended at the top of the output files
    generated at various stages in the pipeline.
    """
    template=textwrap.dedent("""
    ###########################################################
    # File:   %s
    # Author: %s
    # Email:  %s
    # Date:   %s UT
    #
    # Description:
    #
    # %s
    #
    # Command used to generate file:%s
    #
    ###########################################################
    """).lstrip()

    # format description
    desc_raw = " ".join(description.split())
    desc_processed = "\n# ".join(textwrap.wrap(desc_raw, 78))

    # format command
    command_parts = textwrap.wrap(" ".join(sys.argv), 75)
    command_parts = [x.ljust(75) + " \\" for x in command_parts]
    command = "\n# ".join(["\n#"] + command_parts)

    return template % (filename, author, email, datetime.datetime.utcnow(),
                       desc_processed, command)

#--------------------------------------
# Main
#--------------------------------------
args = parse_input()

#--------------------------------------
# Ruffus tasks 
#--------------------------------------
def setup():
    """Create working directories, etc."""
    subdirs = ['01-individual_filtered_reads', '02-combined_filtered_reads']
    for d in [os.path.join('build', x) for x in subdirs]:
        if not os.path.exists(d):
            os.makedirs(d, mode=0755)

@follows(setup)
@transform(args.input_reads, regex(r"^((.*)/)?(.+)\.fastq"),
           r'build/01-individual_filtered_reads/\3.txt', 
           args.spliced_leader, args.min_length)
def parse_reads(infile, outfile, spliced_leader, min_length):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing the feature of interest (SL or polyA tail).
    """
    matches = []

    # limit to matches of size min_length or greater
    s = spliced_leader[-min_length:]

    # compile regex allowing a single mismatch;
    # we will filter based on number of mistmatches again later on a larger
    # portion of the read using the value specified at run-time
    regex = re.compile('|'.join('%s.%s' % (s[:i], s[i+1:])
                           for i in range(len(s))))

    # open fastq file
    fastq = open(infile)
    print("Processing %s" % os.path.basename(infile))

    # filter by length
    for i, read in enumerate(readfq(fastq)):
        seq = read[1][:len(spliced_leader)]

        # check for match
        if re.search(regex, seq) is not None:
            matches.append(seq)

    # save matched reads to file
    with open(outfile, 'w') as fp:
        fp.write("\n".join(matches))

@merge(parse_reads, 
       'build/02-combined_filtered_reads/matching_reads_all_samples.txt')
def filter_reads(input_files, output_file):
    """
    Loads reads from fastq files, filters them, and combined output (for
    now) into a single file.
    """
    # write output
    with open(output_file, 'w') as outfile:
        # add header comment
        header_comment = create_header_comment(
            os.path.basename(output_file),
            """A collection of trimmed and filtered RNA-Seq reads containing at 
             least some portion of the UTR feature of interest in the correct 
             position. This file contains the combined results from all samples
             included in the analysis.""",
            args.author,
            args.email
        )
        outfile.write(header_comment)

        # add header fields
        outfile.write("id, sequence\n")

        for x in input_files:
            with open(x) as infile:
                outfile.write(infile.read() + "\n")

# run pipeline
pipeline_run([filter_reads], verbose=True, multiprocess=12)

