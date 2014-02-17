#!/bin/env python2
# -*- coding: utf-8 -*-
"""
Utranslated Region (UTR) analysis
Keith Hughitt

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

TODO
----
- Generate statistics/plot including:
    - total number of reads (or number of reads mapped to pathogen)
    - number of reads remaining after filtering (reads containing SL)
    - distribution of SL fragment lengths
    - distribution of UTR lengths
    - average number of accepter sites per CDS
- Compare above statistics across samples
- Adding logging (use logging module)

Testing
-------

utr_analysis.py \
    -i "$RAW/tcruzir21/HPGL0258/processed/*.filtered.fastq" \
    -f $REF/tcruzi_clbrener/genome/tc_esmer/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta \
    -g $REF/tcruzi_clbrener/annotation/tc_esmer/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like.gff      \
    -s AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG
"""
import os
import re
import csv
import sys
import glob
import gzip
import pysam
import pandas
import random
import logging
import argparse
import datetime
import StringIO
import textwrap
import subprocess
from ruffus import *
from BCBio import GFF
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

#--------------------------------------
# FASTQ row indices
#--------------------------------------
ID = 0
SEQUENCE = 1
QUALITY = 2

#--------------------------------------
# Non-Ruffus functions
#--------------------------------------
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
        --build-directory tcruzi -m 12 -n 3
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
    parser.add_argument('-d', '--build-directory', required=True,
                        help='Top-level build directory name to use.')
    parser.add_argument('-f', '--fasta-genome', dest='genome', required=True,
                        help='Genome sequence FASTA filepath')
    parser.add_argument('-g', '--gff-annotation', dest='gff', required=True,
                        help='Genome annotation GFF')
    parser.add_argument('-s', '--sl-sequence', dest='spliced_leader', 
                        help='Spliced leader DNA sequence', default=None)
    parser.add_argument('-l', '--left-anchor-reads', dest='anchor_left',
                        help=('Require sequence of interest to be at left '
                              'side of read'), action='store_true')
    parser.add_argument('-r', '--right-anchor-reads', dest='anchor_right',
                        help=('Require sequence of interest to be at right '
                              'side of read'), action='store_true')
    parser.add_argument('-m', '--min-length', default=10, type=int,
                        help='Minimum length of SL match (default=10)')
    parser.add_argument('-n', '--num-mismatches', default=0, type=int,
                        help='Number of mismatches to allow (default=0)')
    parser.add_argument('-w', '--window-size', default=10000, type=int,
                        help=('Number of bases up/downstream of read to look '
                              'for corresponding genes'))
    parser.add_argument('-t', '--num-threads', default=4, type=int,
                        help='Number of threads to use.')
    parser.add_argument('-a', '--author', help='Author contact name', 
                        default='')
    parser.add_argument('-e', '--email', help='Author contact email address',
                        default='')

    # Parse arguments
    args = parser.parse_args()

    # Replace any environmental variables and return args
    args.genome = os.path.expandvars(args.genome)
    args.input_reads = os.path.expandvars(args.input_reads)
    args.gff = os.path.expandvars(args.gff)

    # set defaults for author if none is specified
    if args.author is "":
        args.author = os.getlogin()

    # @TODO Validate input
    return args

def num_lines(filepath):
    """Returns the number of lines in a specified file"""
    if filepath.endswith('.gz'):
        fp = gzip.open(filepath, 'rb')
    else:
        fp = open(filepath)

    # count number of lines
    for i, line in enumerate(fp, 1):
        pass

    fp.close()
    return i

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
        # Modified to preserve full identifier
        # Keith 2014-01-31
        # name, seqs, last = last[1:].partition(" ")[0], [], None
        name, seqs, last = last, [], None
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

def sort_and_index(base_output, num_threads=1):
    """Sorts and indexes .bam files using samtools"""
    sort_cmd = ['samtools', 'sort', '-@', str(num_threads), 
                base_output + ".bam", base_output + "_sorted"]
    loggers[hpgl_id].info(" ".join(sort_cmd))
    subprocess.call(sort_cmd)

    index_cmd = ['samtools', 'index', base_output + '_sorted.bam']
    loggers[hpgl_id].info(" ".join(index_cmd))
    subprocess.call(index_cmd)

#def add_sam_header(filepath, description, author, email, cmd):
    #"""Adds a header to a sam file including some information about how the
    #file was generated, along with some other basic information."""
    #template=textwrap.dedent(""" 
    #@CO File:   %s
    #@CO Author: %s
    #@CO Email:  %s
    #@CO Date:   %s UT
    #@CO Description: %s
    #""").lstrip()

    ## format description
    #desc_raw = " ".join(description.split())
    #desc_processed = "\n@CO ".join(textwrap.wrap(desc_raw, 78))

    ## format command
    #command_parts = textwrap.wrap(cmd, 75)
    #command_parts = [x.ljust(75) + " \\" for x in command_parts]
    #command = "\n@CO ".join(["\n@CO"] + command_parts)

    #return template % (filename, author, email, datetime.datetime.utcnow(),
                       #desc_processed, command)


def run_tophat(output_dir, genome, log_handle, r1, r2="", num_threads=1, 
               max_multihits=20, extra_args=""):
    """
    Uses Tophat to map reads with the specified settings.

    @NOTE 2014/02/06 -- Tophat 2.0.10 fails for some reads when attempting
    to use more than one thread. For now, just perform the mapping
    using a single-thread to be safe...
    """
    # create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o755)

    # build command
    tophat_cmd = (["tophat", "--num-threads", str(num_threads), 
                   "--max-multihits",
                   str(max_multihits)] + extra_args.split() +
                   ["-o", output_dir, genome, r1, r2])

    # run tophat
    log_handle.info(" ".join(tophat_cmd))
    ret = subprocess.call(tophat_cmd)

    # check to see if tophat succeeded
    if ret != 0:
        return ret

    # sort and index bam output using samtools
    sort_and_index(os.path.join(output_dir, 'accepted_hits'), num_threads)
    sort_and_index(os.path.join(output_dir, 'unmapped'), num_threads)

    return 0

def gzip_str(filepath, strbuffer):
    """Takes a StringIO buffer and writes the output to a gzip-compressed
    file."""
    # go to beginning of string buffer
    strbuffer.seek(0)

    # output path
    if filepath.endswith('.gz'):
        outfile = filepath
    else:
        outfile = filepath + '.gz'

    # write contents to a gzip-compressed file
    fp = gzip.open(outfile, 'wb')
    fp.write(strbuffer.read())
    fp.close()

def filter_fastq(infile, outfile, read_ids, log_handle):
    """Takes a filepath to a FASTQ file and returns a new version of the file
    which contains only reads with the specified ids"""
    if infile.endswith('.gz'):
        fastq = gzip.open(infile, 'rb')
    else:
        fastq = open(infile)

    filtered_reads = StringIO.StringIO()

    # get number of reads to be searched
    num_reads = num_lines(infile) / 4
    hundreth = int(round(num_reads / 100))

    # normalize fastq ids and ignore right-part of id row, including the
    # mated pair read number, during comparision
    read_ids = [x.split()[0].strip('@') for x in read_ids]

    # remove redundant entries
    read_ids = list(set(read_ids))

    # read ids may be semi-sorted which can be bad for indexing performance;
    # shuffle to speed things up
    random.shuffle(read_ids)

    log_handle.info(
        "# Filtering fastq for %d matched reads" % len(read_ids)
    )

    # iterate through each entry in fastq file
    for i, read in enumerate(readfq(fastq), 1):
        # log progress
        if i % hundreth == 0:
            log_handle.info('# %2d%%' % round(i / float(hundreth)))

        # normalized entry id
        fastq_id = read[ID].split()[0].strip('@')

        # check to see if it is in the filtered list
        try:
            idx = read_ids.index(fastq_id)
        except ValueError:
            idx = None

        # if it is, add to filtered output
        if idx is not None:
            fastq_entry = [read[ID], read[SEQUENCE], "+", read[QUALITY]]
            filtered_reads.write("\n".join(fastq_entry) + "\n")

            # remove id from list to reduce search space in future iterations
            read_ids.pop(idx)

        # exit loop when all ids have been found
        if len(read_ids) == 0:
            break

    log_handle.info(' # 100%')
    log_handle.info(' # Done filtering')

    # write matching paired-end reads to compressed fastq
    gzip_str(outfile, filtered_reads)
    filtered_reads.close()

#--------------------------------------
# Main
#--------------------------------------

# parse input
args = parse_input()

# analysis type (sl/poly-a)
analysis_type = ('spliced_leader' if args.spliced_leader is not None else
                 'poly-a')

# create a unique build path for specified parameters
subdir = os.path.join(
    args.build_directory,
    analysis_type,
    'mismatches-%d' % args.num_mismatches,
    'minlength-%d' % args.min_length
)
if args.anchor_left:
    subdir = os.path.join(subdir, 'anchor-left')
if args.anchor_right:
    subdir = os.path.join(subdir, 'anchor-right')

# get a list of HPGL ids
input_regex = re.compile(r'.*(HPGL[0-9]+).*')
hpgl_ids = []

for filename in glob.glob(args.input_reads):
    hpgl_id = re.match(input_regex, filename).groups()[0]
    if hpgl_id not in hpgl_ids:
        hpgl_ids.append(hpgl_id)

# create subdirs based on matching parameters
base_dir = os.path.join('build', subdir)

for hpgl_id in hpgl_ids:
    # create build directories
    build_dir = os.path.join('build', subdir, hpgl_id)

    for d in ['fastq', 'ruffus', 'tophat']:
        if not os.path.exists(os.path.join(build_dir, d)):
            os.makedirs(os.path.join(build_dir, d), mode=0o755)

# setup master logger
log_format = '%(asctime)s %(message)s'
date_format = '%Y-%m-%d %I:%M:%S %p'
formatter = logging.Formatter(log_format, datefmt=date_format)

logging.basicConfig(filename=os.path.join(base_dir, 'build.log'),
                    level=logging.INFO,
                    format=log_format,
                    datefmt=date_format)

# log to console as well
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)

logging.info("# Starting UTR Analysis")
logging.info("# Python %s" % sys.version)
logging.info("# Command:\n%s" % " ".join(sys.argv))

# create dictionary of log handlers for sample-specific info
loggers = {}

# setup sample-specific loggers
for hpgl_id in hpgl_ids:
    build_dir = os.path.join('build', subdir, hpgl_id)
    log_file = os.path.join(build_dir, '%s.log' % hpgl_id)
    loggers[hpgl_id] = logging.getLogger(hpgl_id)
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    loggers[hpgl_id].addHandler(handler)


#--------------------------------------
# Ruffus tasks
#--------------------------------------
#
# Input regex explanation:
#
# Ex. "$RAW/tcruzir21/HPGL0121/processed/HPGL0121_R1_filtered.fastq"
#
# \1 - directory
# \2 - HPGLxxxx
# \3 - _anything_between_id_and_read_num_
# \4 - R1/R2
# \5 - _anything_after_read_num_
#
def check_for_bowtie_index():
    """check for bowtie 2 indices and create if needed"""
    genome = os.path.splitext(args.genome)[0]

    # if index exists, stop here
    if os.path.exists('%s.1.bt2' % genome):
        return

    # otherwise, create bowtie2 index
    logging.info("# Building bowtie2 index for %s" % args.genome)
    bowtie_cmd = (['bowtie2-build', args.genome, genome])
    logging.info("# Command:\n" + " ".join(bowtie_cmd))
    ret = subprocess.call(bowtie_cmd)

@follows(check_for_bowtie_index)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq'),
           r'build/%s/\2/ruffus/\2_\4.parse_reads' % subdir,
           r'\2', r'\3', r'\4', r'\5',
           args.spliced_leader, args.min_length)
def parse_reads(input_file, output_file, hpgl_id, file_prefix, read_num, 
                file_suffix, spliced_leader, min_length):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing the feature of interest (SL or polyA tail).

    Output files
    ------------
    There are three possible sets of output files for this function depending
    on whether the input read comes from a mated-pair of reads or a single
    read, and whether (for the case of mated-pair reads) it is the left read
    or right read:

    1. Paired-end, SL suffix in R1
        *_R1_match_R1_with_sl.fastq
        *_R1_match_R1_without_sl.fastq
        *_R1_match_R2.fastq
    2. Paired-end, SL suffix in R2
        *_R2_match_R2_with_sl.fastq
        *_R2_match_R2_without_sl.fastq
        *_R2_match_R1.fastq
    3. Single-end reads
        *_R1_match_with_sl.fastq
        *_R1_match_without_sl.fastq
    """
    # list to keep track of potential SL reads
    matches = []

    # limit to matches of size min_length or greater
    suffix = spliced_leader[-min_length:]

    # Regex position anchors (optional)
    re_prefix = ""
    re_suffix = ""

    if args.anchor_left:
        re_prefix="^"
    if args.anchor_right:
        re_suffix="$"

    # To speed things up, we first filter the reads to find all possible hits
    # by grepping for reads containing at least `min_length` bases of the SL 
    # sequence.

    # If `num_mismatches` is set to 0, only exact matches are
    # allowed. Otherwise a regex is compiled which allows one mismatch at any
    # position in the first `min_length` bases. While this is not ideal (if
    # the user has specified a larger value for `num_mismatches` and more than
    # one occur in this region they will still get filtered out), this
    # speeds things up significantly generally should not result in many real
    # SL hits from being removed.
    if args.num_mismatches == 0:
        read_regex = re.compile(re_prefix + suffix + re_suffix)
    else:
        read_regex = re.compile('|'.join("%s%s.%s%s" % (
            re_prefix, suffix[:i], suffix[i+1:], re_suffix
        ) for i in range(len(suffix))))

    # open fastq file
    fastq = open(input_file)

    # start sample log
    loggers[hpgl_id].info("# Processing %s" % os.path.basename(input_file))

    # total number of reads
    num_reads = num_lines(input_file) / 4
    loggers[hpgl_id].info(
        "# Scanning %d reads for feature of interest (%s: %s)" %
        (num_reads, hpgl_id, read_num)
    )

    # open output string buffer (will write to compressed file later)
    reads_without_sl = StringIO.StringIO()
    reads_with_sl = StringIO.StringIO()

    # Keep track of matched read IDs
    read_ids = []

    # find all reads containing at least `min_length` bases of the feature
    # of interested
    for i, read in enumerate(readfq(fastq)):
        seq = read[SEQUENCE][:len(spliced_leader)]

        # check for match
        match = re.search(read_regex, seq)

        # move on the next read if no match is found
        if match is None:
            continue

        # otherwise add to output fastq
        trimmed_read = [read[ID],
                        read[SEQUENCE][match.end():],
                        "+",
                        read[QUALITY][match.end():]]
        reads_without_sl.write("\n".join(trimmed_read) + "\n")

        # also save untrimmed read (for finding false SL hits)
        untrimmed_read = [read[ID],
                          read[SEQUENCE],
                          "+",
                          read[QUALITY]]
        reads_with_sl.write("\n".join(untrimmed_read) + "\n")

        # save id
        read_ids.append(read[ID])

    # log numbers
    loggers[hpgl_id].info(
        "# Found %d reads with possible feature of interest (%s)" %
        (len(read_ids), read_num)
    )

    # Paired-end reads
    if os.path.isfile(input_file.replace('_R1_', '_R2_')):
        paired_end = True

        # matching reads
        output_base = 'build/%s/%s/fastq/possible_sl_reads/%s_%s_match' % (
            subdir, hpgl_id, hpgl_id, read_num 
        )
        output_with_sl = "%s_%s_with_sl.fastq" % (output_base, read_num)
        output_without_sl = "%s_%s_without_sl.fastq" % (output_base, read_num)

        # mated reads
        read_num_other = "R1" if read_num == "R2" else "R2"
        input_mated_reads = input_file.replace(read_num, read_num_other)
        output_mated_reads = "%s_%s.fastq" % (output_base, read_num_other)

    # Single-end reads
    else:
        paired_end = False

        output_base = 'build/%s/%s/fastq/possible_sl_reads/%s_R1' % (
            subdir, hpgl_id, hpgl_id
        )
        output_with_sl = "%s_with_sl.fastq" % output_base
        output_without_sl = "%s_without_sl.fastq" % output_base

    # Create output directory
    output_dir = os.path.dirname(output_base)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o755)

    # Save complete (untrimmed) reads containing a portion of the SL sequence
    # as well. By mapping these reads to the genome we can find false hits -- 
    # i.e. reads that contain a part of the SL sequence that is not actually 
    # from the SL.

    # write trimmed and untrimmed reads to fastq.gz
    gzip_str(output_without_sl, reads_without_sl)
    gzip_str(output_with_sl, reads_with_sl)

    # clean up
    fastq.close()
    reads_without_sl.close()
    reads_with_sl.close()

    loggers[hpgl_id].info(("# Finished processing %s" %
                           os.path.basename(input_file)))

    # For PE reads, grab reads from matching pair
    if paired_end:
        loggers[hpgl_id].info(("# Retrieving mated pair reads for %s" % 
                               os.path.basename(input_mated_reads)))
        filter_fastq(input_mated_reads, output_mated_reads, read_ids,
                     loggers[hpgl_id])

    # Let Ruffus know we are done
    open(output_file, 'w').close()

@transform(parse_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).parse_reads'),
           r'\1/\2_\3.remove_false_hits',
           r'\2', r'\3')
def remove_false_hits(input_file, output_file, hpgl_id, read_num):
    output_dir = ('build/%s/%s/tophat/%s_false_hits' % (subdir, hpgl_id,
                                                        read_num))
    genome = os.path.splitext(args.genome)[0]

    # input read base directory
    basedir = 'build/%s/%s/fastq' % (subdir, hpgl_id)

    # R1 filepath (including matched SL sequence)
    if (read_num == 'R1'):
        r1_filepath = ('%s/possible_sl_reads/%s_R1_match_R1_with_sl.fastq.gz' %
                       (basedir, hpgl_id))

        # R2 filepath (for PE reads)
        r2_filepath = r1_filepath.replace('R1_with_sl', 'R2')

        # If SE, set filepath to empty string
        if not os.path.exists(r2_filepath):
            r2_filepath = ""
    # R2
    else:
        r2_filepath = ('%s/possible_sl_reads/%s_R2_match_R2_with_sl.fastq.gz' %
                       (basedir, hpgl_id))
        r1_filepath = r2_filepath.replace('R2_with_sl', 'R1')

    # Map reads using Tophat
    # @TODO parameterize extra_args (except for --no-mixed in this case) to
    # allow for easier customization
    loggers[hpgl_id].info(("# Mapping full reads containing feature of "
                           "interest to find false hits (reads that correspond"
                           " to actual features in the genome"))
    ret = run_tophat(output_dir, genome, loggers[hpgl_id],
                     r1_filepath, r2_filepath,
                     extra_args='--mate-inner-dist 170 --no-mixed')

    # Make sure tophat succeeded
    if ret != 0:
        logging.error(
            "# Error running tophat 1/2! %s (%s)" % (hpgl_id, read_num)
        )
        sys.exit()

    # Get ids of actual SL-containing reads (those that failed to map when the
    # SL sequence was included).
    sam = pysam.Samfile(os.path.join(output_dir, 'unmapped.bam'), 'rb')
    good_ids = [x.qname for x in sam]

    # number of reads before filtering
    num_reads_before = num_lines(r1_filepath) / 4
    loggers[hpgl_id].info(("# Removing %d false hits (%d total)" % (
        num_reads_before - good_ids, num_reads_before
    )))

    # Create true hits directory
    hits_dir = os.path.join(basedir, 'actual_sl_reads')
    if not os.path.exists(hits_dir):
        os.makedirs(hits_dir, mode=0o755)

    # Create filtered versions of R1 (and R2) fastq files with only the un-
    # mapped reads

    # Filter R1 reads
    r1_infile = r1_filepath.replace('with_sl', 'without_sl')
    r1_outfile = r1_infile.replace('possible', 'actual')
    loggers[hpgl_id].info(
        "# Filtering matched reads to remove false hits (R1)"
    )
    filter_fastq(r1_infile, r1_outfile, good_ids, loggers[hpgl_id])

    # Filter R2 reads
    if r2_filepath != "":
        r2_infile = r2_filepath.replace('with_sl', 'without_sl')
        r2_outfile = r2_infile.replace('possible', 'actual')
        loggers[hpgl_id].info(
            "# Filtering matched reads to remove false hits (R2)"
        )
        filter_fastq(r2_infile, r2_outfile, good_ids, loggers[hpgl_id])

    loggers[hpgl_id].info("# Finished removing false hits.")

    # Let Ruffus know we are done
    open(output_file, 'w').close()

@transform(remove_false_hits,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).remove_false_hits'),
           r'\1/\2_\3.map_sl_reads',
           r'\2', r'\3')
def map_sl_reads(input_file, output_file, hpgl_id, read_num):
    """Maps the filtered spliced-leader containing reads back to the genome"""
    output_dir = 'build/%s/%s/tophat/%s_sl_reads' % (subdir, hpgl_id, read_num)
    genome = os.path.splitext(args.genome)[0]

    loggers[hpgl_id].info("# Mapping filtered reads back to genome")

    # input read base directory
    basedir = 'build/%s/%s/fastq' % (subdir, hpgl_id)

    # R1 filepath (including matched SL sequence)
    if (read_num == 'R1'):
        r1_filepath = (
            '%s/actual_sl_reads/%s_R1_match_R1_without_sl.fastq.gz' %
            (basedir, hpgl_id)
        )

        # R2 filepath (for PE reads)
        r2_filepath = r1_filepath.replace('R1_with_sl', 'R2')

        # If SE, set filepath to empty string
        if not os.path.exists(r2_filepath):
            r2_filepath = ""
    # R2
    else:
        r2_filepath = (
            '%s/actual_sl_reads/%s_R2_match_R2_without_sl.fastq.gz' %
            (basedir, hpgl_id)
        )
        r1_filepath = r2_filepath.replace('R2_with_sl', 'R1')

    # Map reads using Tophat
    #  --no-mixed ?
    ret = run_tophat(output_dir, genome, loggers[hpgl_id],
                 r1_filepath, r2_filepath,
                 extra_args='--mate-inner-dist 170 --transcriptome-max-hits 1')

    # Make sure tophat succeeded
    if ret != 0:
        logging.error(
            "# Error running tophat 2/2! %s (%s)" % (hpgl_id, read_num)
        )
        sys.exit()

    loggers[hpgl_id].info("# Finished mapping hits to genome")

    # Let Ruffus know we are done
    open(output_file, 'w').close()

@collate(map_sl_reads,
       regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).map_sl_reads'),
       r'\1/\2_.compute_coordinates')
def compute_coordinates(input_files, output_file):
    """Maps the filtered spliced-leader containing reads back to the genome.

    References
    ----------
    * http://www.cgat.org/~andreas/documentation/pysam/api.html#pysam.Samfile
    * http://biopython.org/wiki/GFF_Parsing
    * http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html
    """
    # load GFF
    gff_fp = open(args.gff)

    logging.info("# Computing coordinates for mapped hits")

    # Get chromosomes from GFF file
    # @TODO: Generalize for other species
    chromosomes = {}
    for rec in GFF.parse(gff_fp):
        #if rec.id.startswith('TcChr'):
        if len(rec.features) > 0 and rec.features[0].type == 'chromosome':
            chromosomes[rec.id] = rec

    # Create a dictionary to keep track of the splice acceptor site
    # locations and frequencies
    # A nested dictionary will be used such that a given SL site can be
    # indexed using the syntax:
    # results[chr_num][gene_id][acceptor_site_offset]
    results = {}

    # Bam inputs
    input_globstr = ('build/%s/*/tophat/*_sl_reads/accepted_hits_sorted.bam' %
                     subdir)

    # TESTING 2014/02/13
    # Let's see what is not matching
    not_matched = 'build/%s/no_matches.csv' % (subdir)
    debug_writer = csv.writer(open(not_matched, 'w'))
    debug_writer.writerow(['read_id', 'chromosome', 'strand', 'position'])

    # Itereate over mapped reads
    for filepath in glob.glob(input_globstr):
        # get hpgl id
        hpgl_id = re.match('.*(HPGL[0-9]+).*', x).groups()[0]

        # file to save results for individual sample
        base_dir = filepath[:re.search('tophat', filepath).start()]
        sample_hits = os.path.join(base_dir, hpgl_id + '_sl_coordinates.csv')
        sample_csv_writer = csv.writer(open(sample_hits, 'w'))
        sample_csv_writer.writerow(['read_id', 'gene_id', 'chromosome', 
                                    'strand', 'position', 'distance'])

        # keep track of how many reads were found in the expected location
        num_good = 0
        num_bad = 0

        # open sam file
        sam = pysam.Samfile(filepath, 'rb')

        # Keep track of read id so we only count each one once
        read_ids = []

        # Get coordinate and strand for each read in bam file
        for read in sam:
            # if we have already counted the read, stop here
            if read.qname in read_ids:
                continue
            read_ids.append(read.qname)

            pos = read.pos
            chromosome = sam.getrname(read.tid)
            strand = -1 if read.is_reverse else 1

            # Find nearest gene
            # 1. Get genes within +/- N bases of location (if any)
            # 2. Find closest match
            if strand == 1:
                # For positive-strand sites, search region just downstream 
                # of splice-site for genes
                subseq = chromosomes[chromosome][pos:pos + args.window_size]
                gene_start = 'start'
                offset = 0
            else:
                # For negative-strand sites, search region just upstream
                # of splice-site for genes
                subseq = chromosomes[chromosome][pos - args.window_size:pos]
                gene_start = 'end'
                offset = args.window_size

            # If there are no nearby genes, stop here
            if len(subseq.features) == 0:
                debug_writer.writerow([read.qname, chromosome, strand, pos])
                num_bad = num_bad + 1
                continue

            num_good = num_good + 1

            # Otherwise find closest gene to the acceptor site
            closest_gene = None
            closest_dist = float('inf')

            for f in subseq.features:
                dist = abs(offset - int(getattr(f.location, gene_start)))
                if dist < closest_dist:
                    closest_gene = f.id
                    closest_dist = dist

            # Add to output dictionary
            if not chromosome in results:
                results[chromosome] = {}
            if not closest_gene in results[chromosome]:
                results[chromosome][closest_gene] = {}

            # Add entry to sample output csv
            sample_csv_writer.writerow([
                read.qname, closest_gene, chromosome, strand, pos, closest_dist
            ])

            # Increment SL site count and save distance from gene
            if not pos in results[chromosome][closest_gene]:
                results[chromosome][closest_gene][pos] = {
                    "count": 1,
                    "dist": closest_dist
                }
            else:
                results[chromosome][closest_gene][pos]['count'] += 1

        # record number of good and bad reads
        loggers[hpgl_id].info(
            "# Found %d reads with feature in expected location" % num_good)
        loggers[hpgl_id].info(
            "# Found %d reads with feature in incorrect location" % num_bad)

    # Output filepath
    coordinates_output = 'build/%s/sl_coordinates.csv' % (subdir)
    fp = open(coordinates_output, 'w')

    # Write csv header
    header = create_header_comment(os.path.basename(coordinates_output),
                                   "Spliced leader acceptor site coordinates",
                                   args.author, args.email)
    fp.write(header)

    # Write header to output
    writer = csv.writer(fp)
    writer.writerow(['gene', 'chromosome', 'location', 'distance', 
                     'count'])

    # write output to csv
    for chrnum in results:
        for gene_id in results[chrnum]:
            for acceptor_site in results[chrnum][gene_id]:
                writer.writerow([
                    gene_id, chrnum, acceptor_site,
                    results[chrnum][gene_id][acceptor_site]['dist'],
                    results[chrnum][gene_id][acceptor_site]['count']
                ])

    logging.info("# Finished!")

    # clean up
    gff_fp.close()
    fp.close()

#@follows(filter_reads)
#def compute_read_statistics():
    #"""Computes some basic stats about the distribution of the UTR feature
    #in the RNA-Seq reads"""
    #subdir = 'mismatches-%d_minlength-%d' % (args.num_mismatches, args.min_length)
    #fp = open('build/02-combined_filtered_reads/%s/matching_reads_all_samples.csv' %
              #subdir)

    ## skip comments
    #line = fp.readline()
    #while line.startswith('#'):
        #line = fp.readline()
        #continue

    ## col names
    #colnames = line.strip().split(',')

    ## load csv into a pandas dataframe
    #df = pandas.read_csv(fp, header=None, names=colnames)
    #fp.close()

    ## summary statistics
    #print("SL read length distribution:")
    #print(df.groupby('end').size())
    #print(df['end'].describe())

    ## plot a histogram of the SL lengths captured in the RNA-Seq reads
    #df.hist(column='end', bins=len(args.spliced_leader) - args.min_length)
    #title_text = "SL fragment length distribution ($N=%d$, $\mu=%f$)"
    #plt.title(title_text % (df['end'].size, df['end'].mean()))
    ##plt.savefig('output/figures/%s/sl_length_distribution.png' % subdir)
    #plt.savefig('sl_length_distribution.png' % subdir)

# run pipeline
if __name__ == "__main__":
    pipeline_run([compute_coordinates], logger=logging.getLogger(''),
                 multiprocess=args.num_threads)
    pipeline_printout_graph("utr_analysis_flowchart.png", "png",
                            [compute_coordinates])

