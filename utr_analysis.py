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
- For reads where feature was found in R2, discard R1 and just map R2
  (see RIP notes)
- Generate statistics/plot including:
    - total number of reads (or number of reads mapped to pathogen)
    - number of reads remaining after filtering (reads containing SL)
    - distribution of SL fragment lengths
    - distribution of UTR lengths
    - average number of accepter sites per CDS
- Compare above statistics across samples

Testing
-------

utr_analysis.py \
    -i "$RAW/tcruzir21/HPGL0258/processed/*.filtered.fastq" \
    -f $REF/tcruzi_clbrener/genome/tc_esmer/TriTrypDB-7.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta \
    -g $REF/tcruzi_clbrener/annotation/tc_esmer/TriTrypDB-7.0_TcruziCLBrenerEsmeraldo-like.gff      \
    -s AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG
"""
import os
import re
import csv
import sys
import glob
import gzip
import pysam
import random
import logging
import argparse
import datetime
import StringIO
import textwrap
import subprocess
from ruffus import *
from BCBio import GFF

#--------------------------------------
# FASTQ row indices
#--------------------------------------
ID_IDX = 0
SEQUENCE_IDX = 1
QUALITY_IDX = 2

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
        -f TriTrypDB-7.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta \\
        -g TrypDB-7.0_TcruziCLBrenerEsmeraldo-like.gff             \\
        --build-directory build/tcruzi --min-sl-length 12
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
                        help='Directory to save output to')
    parser.add_argument('-f', '--fasta-genome', dest='genome', required=True,
                        help='Genome sequence FASTA filepath')
    parser.add_argument('-g', '--gff-annotation', dest='gff', required=True,
                        help='Genome annotation GFF')
    parser.add_argument('-s', '--sl-sequence', dest='spliced_leader',
                        required=True, help='Spliced leader DNA sequence', 
                        default=None)
    parser.add_argument('-e', '--exclude-internal-matches',
                        help=('Only allow matches with the feature at the '
                              'expected end of a read (upstream for SL and'
                              'downstream for Poly-A tail)'),
                        action='store_true')
    parser.add_argument('-m', '--min-sl-length', default=10, type=int,
                        help='Minimum length of SL match (default=10)')
    parser.add_argument('-p', '--min-polya-length', default=10, type=int,
                        help='Minimum length of Poly-A match (default=10)')
    parser.add_argument('-w', '--window-size', default=10000, type=int,
                        help=('Number of bases up/downstream of read to look '
                              'for corresponding genes'))
    parser.add_argument('-t', '--num-threads', default=4, type=int,
                        help='Number of threads to use.')
    parser.add_argument('-a', '--author', help='Author contact name', 
                        default='')
    parser.add_argument('-c', '--email', help='Author contact email address',
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

def run_command(cmd, log_handle):
    """Runs a command and logs the output to a specified log handle"""
    log_handle.info(cmd)

    process = subprocess.Popen(cmd.split(" "),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if stdout:
        log_handle.info(stdout)
    if stderr:
        log_handle.error(stderr)

    return process.returncode

def sort_and_index(base_output, log_handle):
    """Sorts and indexes .bam files using samtools"""
    # sort bam
    sort_cmd = 'samtools sort %s %s' % (
        base_output + ".bam", base_output + "_sorted")
    run_command(sort_cmd, log_handle)

    # index bam
    index_cmd = 'samtools index %s' % (base_output + '_sorted.bam')
    run_command(index_cmd, log_handle)

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
    cmd = "tophat --num-threads %d --max-multihits %d %s -o %s %s %s %s" % (
           num_threads, max_multihits, extra_args, output_dir, genome, r1, r2)

    # run tophat
    ret = run_command(cmd, log_handle)

    # check to see if tophat succeeded
    if ret != 0:
        return ret

    # sort and index bam output using samtools
    sort_and_index(os.path.join(output_dir, 'accepted_hits'), log_handle)
    sort_and_index(os.path.join(output_dir, 'unmapped'), log_handle)

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

def filter_fastq(infile1, infile2, outfile1, outfile2, read_ids, log_handle):
    """Takes a filepath to a FASTQ file and returns a new version of the file
    which contains only reads with the specified ids"""
    if infile1.endswith('.gz'):
        fastq1 = gzip.open(infile1, 'rb')
        fastq2 = gzip.open(infile2, 'rb')
    else:
        fastq1 = open(infile1)
        fastq2 = open(infile2)

    filtered_reads1 = StringIO.StringIO()
    filtered_reads2 = StringIO.StringIO()

    # get number of reads to be searched
    num_reads = num_lines(infile1) / 4
    hundreth = int(round(num_reads / 100))
    show_progress = len(read_ids) >= 10000

    # normalize fastq ids and ignore right-part of id row, including the
    # mated pair read number, during comparision
    read_ids = [x.split()[0].strip('@') for x in read_ids]

    # remove redundant entries
    read_ids = list(set(read_ids))

    # read ids may be semi-sorted which can be bad for indexing performance;
    # shuffle to speed things up
    random.shuffle(read_ids)

    log_handle.info("# Filtering fastq for %d matched reads" % len(read_ids))

    # mated reads handle
    mated_reads = readfq(fastq2)

    # iterate through each entry in fastq file
    for i, read in enumerate(readfq(fastq1), 1):
        mated_read = mated_reads.next()

        # log progress for large numbers of records
        if show_progress and (i % hundreth) == 0:
            log_handle.info('# %2d%%' % round(i / float(hundreth)))

        # normalized entry id
        fastq_id = read[ID_IDX].split()[0].strip('@')

        # check to see if it is in the filtered list
        try:
            idx = read_ids.index(fastq_id)
        except ValueError:
            idx = None

        # if it is, add to filtered output
        if idx is not None:
            # read 1
            fastq_entry1 = [read[ID_IDX], read[SEQUENCE_IDX], "+", 
                            read[QUALITY_IDX]]
            filtered_reads1.write("\n".join(fastq_entry1) + "\n")

            # read 2
            fastq_entry2 = [mated_read[ID_IDX], mated_read[SEQUENCE_IDX], "+", 
                            mated_read[QUALITY_IDX]]
            filtered_reads2.write("\n".join(fastq_entry2) + "\n")


            # remove id from list to reduce search space in future iterations
            read_ids.pop(idx)

        # exit loop when all ids have been found
        if len(read_ids) == 0:
            break

    log_handle.info(' # 100%')
    log_handle.info(' # Done filtering')

    # write matching paired-end reads to compressed fastq
    gzip_str(outfile1, filtered_reads1)
    gzip_str(outfile2, filtered_reads2)
    filtered_reads1.close()
    filtered_reads2.close()

def get_next_log_name(base_name):
    """Returns a filepath for the next highest log number"""
    if not os.path.exists(base_name):
        return base_name
    else:
        log_nums = [int(x.split('.').pop()) 
                        for x in glob.glob("%s.*" % base_name)]
        next_log_num = max([0] + log_nums) + 1
        return "%s.%d" % (base_name, next_log_num)

#--------------------------------------
# Main
#--------------------------------------

# parse input
args = parse_input()

# analysis type (sl/poly-a)
#analysis_type = ('spliced_leader' if args.spliced_leader is not None else
#                 'poly-a')

# create a unique build path for specified parameters

# SL sub-directory
sl_build_dir = os.path.join(
    args.build_directory,
    'spliced_leader',
    'minlength-%d' % args.min_sl_length,
    'anchored' if args.exclude_internal_matches else 'unanchored'
)

# poly-A tail sub-directory
polya_build_dir = os.path.join(
    args.build_directory,
    'poly-a',
    'minlength-%d' % args.min_polya_length,
    'anchored' if args.exclude_internal_matches else 'unanchored'
)

# get a list of HPGL ids
input_regex = re.compile(r'.*(HPGL[0-9]+).*')
hpgl_ids = []

# currently processing
for filename in glob.glob(args.input_reads):
    hpgl_id = re.match(input_regex, filename).groups()[0]
    if hpgl_id not in hpgl_ids:
        hpgl_ids.append(hpgl_id)

# list of ids including previously processed samples (used for final step)
hpgl_ids_all = hpgl_ids
input_globstr = (
    '%s/*/tophat/*_sl_reads/accepted_hits_sorted.bam' % sl_build_dir
)
for filepath in glob.glob(input_globstr):
    # get hpgl id
    hpgl_ids_all.append(re.match('.*(HPGL[0-9]+).*', filepath).groups()[0])

# create subdirs based on matching parameters
for hpgl_id in hpgl_ids:
    # create build directories
    for base_dir in [sl_build_dir, polya_build_dir]:
        for sub_dir in ['fastq', 'ruffus', 'tophat']:
            outdir = os.path.join(base_dir, hpgl_id, sub_dir)
            if not os.path.exists(outdir):
                os.makedirs(outdir, mode=0o755)

# setup master logger
log_format = '%(asctime)s %(message)s'
date_format = '%Y-%m-%d %I:%M:%S %p'
formatter = logging.Formatter(log_format, datefmt=date_format)

# determine log name to use
master_log = get_next_log_name(os.path.join(args.build_directory, 'build.log'))

logging.basicConfig(filename=master_log,
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
for hpgl_id in hpgl_ids_all:
    loggers[hpgl_id] = {}

    for analysis in ['SL', 'PolyA']:
        loggers[hpgl_id][analysis] = {}

        for read_num in ['R1', 'R2']:
            build_dir = sl_build_dir if analysis == 'SL' else polya_build_dir

            sample_log_name = get_next_log_name(
                os.path.join(build_dir, hpgl_id, '%s_%s_%s.log' % (hpgl_id,
                    analysis, read_num))
            )
            loggers[hpgl_id][analysis][read_num] = logging.getLogger(
                hpgl_id + analysis + read_num
            )
            handler = logging.FileHandler(sample_log_name)
            handler.setFormatter(formatter)
            loggers[hpgl_id][analysis][read_num].addHandler(handler)

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

def check_for_genome_fasta():
    """Checks to make sure the genome fasta exists and is available as a .fa
    file; required by Tophat."""
    genome = os.path.splitext(args.genome)[0]

    # if index exists, stop here
    if os.path.exists('%s.fa' % genome):
        return
    elif os.path.exists('%s.fasta' % genome):
        # if .fasta file exists, but not .fa, create a symlink
        logging.info("# Creating symlink to %s for Tophat" % args.genome)
        os.symlink(genome + '.fasta', genome + '.fa')
    else:
        raise IOError("Missing genome file")

@follows(check_for_bowtie_index)
@follows(check_for_genome_fasta)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq'),
           r'%s/\2/ruffus/\2_\4.find_sl_reads' % sl_build_dir,
           r'\2', r'\3', r'\4', r'\5',
           args.spliced_leader, args.min_sl_length)
def find_sl_reads(input_file, output_file, hpgl_id, file_prefix, read_num, 
                file_suffix, spliced_leader, min_sl_length):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing a subsequence of the spliced leader.

    Output files
    ------------
    There are three possible sets of output files for this function depending
    on whether the input read comes from a mated-pair of reads or a single
    read, and whether (for the case of mated-pair reads) it is the left read
    or right read:

    1. SL suffix in R1
        *_R1_1_with_sl.fastq
        *_R1_1_without_sl.fastq
        *_R1_2.fastq
    2. SL suffix in R2
        *_R2_2_with_sl.fastq
        *_R2_2_without_sl.fastq
        *_R2_1.fastq
    """
    # sample log
    log = loggers[hpgl_id]['SL'][read_num]

    # list to keep track of potential SL reads
    matches = []

    # output string buffers and filepaths
    output_base = '%s/%s/fastq/possible_sl_reads/%s_%s' % (
        sl_build_dir, hpgl_id, hpgl_id, read_num 
    )
    output_with_sl = "%s_%s_with_sl.fastq" % (output_base, read_num[-1])
    output_without_sl = "%s_%s_without_sl.fastq" % (output_base, read_num[-1])

    # mated reads
    read_num_other = "R1" if read_num == "R2" else "R2"
    input_file_mated = input_file.replace(read_num, read_num_other)
    output_mated_reads = "%s_%s.fastq" % (output_base, read_num_other[-1])

    # Determine strings to match in reads
    if args.exclude_internal_matches:
        read_regex_str = '|'.join(["^" + spliced_leader[-x:] for x in
            range(min_sl_length, len(spliced_leader) + 1)])
    else:
        read_regex_str = spliced_leader[-min_sl_length:]

    # Compile regular expression
    read_regex = re.compile(read_regex_str)

    # total number of reads
    num_reads = num_lines(input_file) / 4

    # Start sample log
    log.info("# Processing %s" % os.path.basename(input_file))
    log.info("# Scanning %d reads for spliced leader sequence" % (num_reads))
    log.info("# Using Regex patten:\n %s" % read_regex_str)

    # open output string buffer (will write to compressed file later)
    reads_without_sl = StringIO.StringIO()
    reads_with_sl = StringIO.StringIO()
    mated_reads_buffer = StringIO.StringIO()

    # Keep track of matched read IDs
    read_ids = []

    # find all reads containing at least `min_sl_length` bases of the feature
    # of interested
    fastq = open(input_file)
    fastq_mated = open(input_file_mated)

    # iterate over mated reads at same time
    mated_reads = readfq(fastq_mated)

    for i, read in enumerate(readfq(fastq)):
        # get mated read
        mated_read = mated_reads.next()

        # check for match
        match = re.search(read_regex, read[SEQUENCE_IDX])

        # move on the next read if no match is found
        if match is None:
            continue

        # otherwise add to output fastq
        trimmed_read = [read[ID_IDX],
                        read[SEQUENCE_IDX][match.end():],
                        "+",
                        read[QUALITY_IDX][match.end():]]
        reads_without_sl.write("\n".join(trimmed_read) + "\n")

        # Also save complete (untrimmed) reads containing a portion of the SL 
        # sequence as well. By mapping these reads to the genome we can find 
        # false hits; i.e. reads that contain a part of the SL sequence that 
        # is not actually from the SL.
        untrimmed_read = [read[ID_IDX],
                          read[SEQUENCE_IDX],
                          "+",
                          read[QUALITY_IDX]]
        reads_with_sl.write("\n".join(untrimmed_read) + "\n")

        # paired-end reads
        untrimmed_mated_read = [mated_read[ID_IDX],
                                mated_read[SEQUENCE_IDX],
                                "+",
                                mated_read[QUALITY_IDX]]
        mated_reads_buffer.write("\n".join(untrimmed_mated_read) + "\n")

        # save id
        read_ids.append(read[ID_IDX])

    # log numbers
    log.info("# Found %d reads with possible spliced leader fragment" % 
             len(read_ids))

    # Create output directory
    output_dir = os.path.dirname(output_base)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o755)

    # write trimmed and untrimmed reads to fastq.gz
    gzip_str(output_without_sl, reads_without_sl)
    gzip_str(output_with_sl, reads_with_sl)
    gzip_str(output_mated_reads, mated_reads_buffer)

    # clean up
    fastq.close()
    fastq_mated.close()
    reads_without_sl.close()
    reads_with_sl.close()
    mated_reads_buffer.close()

    log.info("# Finished processing %s" % os.path.basename(input_file))

    # Let Ruffus know we are done
    open(output_file, 'w').close()

@follows(check_for_bowtie_index)
@follows(check_for_genome_fasta)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq'),
           r'%s/\2/ruffus/\2_\4.find_polya_reads' % polya_build_dir,
           r'\2', r'\3', r'\4', r'\5',
           args.min_polya_length)
def find_polya_reads(input_file, output_file, hpgl_id, file_prefix, read_num, 
                    file_suffix, min_polya_length):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing the a putative poly-A tail sequence.

    @TODO: look for T's on left, or A's on right

    Output files
    ------------
    There are three possible sets of output files for this function depending
    on whether the input read comes from a mated-pair of reads or a single
    read, and whether (for the case of mated-pair reads) it is the left read
    or right read:

    1. SL suffix in R1
        *_R1_1_with_polya.fastq
        *_R1_1_without_polya.fastq
        *_R1_2.fastq
    2. SL suffix in R2
        *_R2_2_with_pola.fastq
        *_R2_2_without_polya.fastq
        *_R2_1.fastq
    """
    # sample log
    log = loggers[hpgl_id]['PolyA'][read_num]

    # list to keep track of potential SL reads
    matches = []

    # output string buffers and filepaths
    output_base = '%s/%s/fastq/possible_polya_reads/%s_%s' % (
        polya_build_dir, hpgl_id, hpgl_id, read_num 
    )
    output_with_polya = "%s_%s_with_polya.fastq" % (output_base, read_num[-1])
    output_without_polya = "%s_%s_without_polya.fastq" % (output_base, read_num[-1])

    # mated reads
    read_num_other = "R1" if read_num == "R2" else "R2"
    input_file_mated = input_file.replace(read_num, read_num_other)
    output_mated_reads = "%s_%s.fastq" % (output_base, read_num_other[-1])

    # Compile regular expression
    read_regex_str = 'T{%d,}|A{%d,}' % (min_polya_length, min_polya_length)

    if args.exclude_internal_matches:
        read_regex_str = "^" + read_regex_str + "$"

    # total number of reads
    num_reads = num_lines(input_file) / 4

    # Start sample log
    log.info("# Processing %s" % os.path.basename(input_file))
    log.info("# Scanning %d reads for Poly-A tail" % (num_reads))
    log.info("# Using Regex patten:\n %s" % read_regex_str)

    # open output string buffer (will write to compressed file later)
    reads_without_polya = StringIO.StringIO()
    reads_with_polya = StringIO.StringIO()
    mated_reads_buffer = StringIO.StringIO()

    # Keep track of matched read IDs
    read_ids = []

    # find all reads containing at least `min_polya_length` bases of the feature
    # of interested
    fastq = open(input_file)
    fastq_mated = open(input_file_mated)

    # iterate over mated reads at same time
    mated_reads = readfq(fastq_mated)

    for i, read in enumerate(readfq(fastq)):
        # get mated read
        mated_read = mated_reads.next()

        # check for match
        match = re.search(read_regex, read[SEQUENCE_IDX])

        # move on the next read if no match is found
        if match is None:
            continue

        # otherwise add to output fastq
        trimmed_read = [read[ID_IDX],
                        read[SEQUENCE_IDX][match.end():],
                        "+",
                        read[QUALITY_IDX][match.end():]]
        reads_without_polya.write("\n".join(trimmed_read) + "\n")

        # Also save complete (untrimmed) reads containing a portion of the SL 
        # sequence as well. By mapping these reads to the genome we can find 
        # false hits; i.e. reads that contain a part of the SL sequence that 
        # is not actually from the SL.
        untrimmed_read = [read[ID_IDX],
                          read[SEQUENCE_IDX],
                          "+",
                          read[QUALITY_IDX]]
        reads_with_polya.write("\n".join(untrimmed_read) + "\n")

        # paired-end reads
        untrimmed_mated_read = [mated_read[ID_IDX],
                                mated_read[SEQUENCE_IDX],
                                "+",
                                mated_read[QUALITY_IDX]]
        mated_reads_buffer.write("\n".join(untrimmed_mated_read) + "\n")

        # save id
        read_ids.append(read[ID_IDX])

    # log numbers
    log.info("# Found %d reads with possible poly-A tail fragment" % 
             len(read_ids))

    # Create output directory
    output_dir = os.path.dirname(output_base)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o755)

    # write trimmed and untrimmed reads to fastq.gz
    gzip_str(output_without_polya, reads_without_polya)
    gzip_str(output_with_polya, reads_with_polya)
    gzip_str(output_mated_reads, mated_reads_buffer)

    # clean up
    fastq.close()
    fastq_mated.close()
    reads_without_polya.close()
    reads_with_polya.close()
    mated_reads_buffer.close()

    log.info("# Finished processing %s" % os.path.basename(input_file))

    # Let Ruffus know we are done
    open(output_file, 'w').close()

@transform(find_sl_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).find_sl_reads'),
           r'\1/\2_\3.remove_false_hits',
           r'\2', r'\3')
def remove_false_hits(input_file, output_file, hpgl_id, read_num):
    output_dir = ('%s/%s/tophat/%s_false_hits' % (sl_build_dir, hpgl_id,
                                                        read_num))
    genome = os.path.splitext(args.genome)[0]

    # input read base directory
    basedir = '%s/%s/fastq' % (sl_build_dir, hpgl_id)

    # R1 filepath (including matched SL sequence)
    if (read_num == 'R1'):
        r1_filepath = ('%s/possible_sl_reads/%s_R1_1_with_sl.fastq.gz' %
                       (basedir, hpgl_id))
        r2_filepath = r1_filepath.replace('1_with_sl', '2')
    # R2
    else:
        r2_filepath = ('%s/possible_sl_reads/%s_R2_2_with_sl.fastq.gz' %
                       (basedir, hpgl_id))
        r1_filepath = r2_filepath.replace('2_with_sl', '1')

    # Map reads using Tophat
    # @TODO parameterize extra_args (except for --no-mixed in this case) to
    # allow for easier customization
    loggers[hpgl_id]['SL'][read_num].info(
        "# Mapping full reads containing feature of interest to find false\n"
        "# hits (reads that correspond to actual features in the genome"
    )
    ret = run_tophat(output_dir, genome, loggers[hpgl_id]['SL'][read_num],
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
    loggers[hpgl_id]['SL'][read_num].info(
        "# Removing %d false hits (%d total)" %
        (num_reads_before - len(good_ids), num_reads_before)
    )

    # Create true hits directory
    hits_dir = os.path.join(basedir, 'actual_sl_reads')
    if not os.path.exists(hits_dir):
        os.makedirs(hits_dir, mode=0o755)

    # Create filtered versions of R1 (and R2) fastq files with only the un-
    # mapped reads
    loggers[hpgl_id]['SL'][read_num].info(
        "# Filtering matched reads to remove false hits"
    )

    # Filepaths
    r1_infile = r1_filepath.replace('with_sl', 'without_sl')
    r2_infile = r2_filepath.replace('with_sl', 'without_sl')
    r1_outfile = r1_infile.replace('possible', 'actual')
    r2_outfile = r2_infile.replace('possible', 'actual')

    filter_fastq(r1_infile, r2_infile, r1_outfile, r2_outfile,
                 good_ids, loggers[hpgl_id]['SL'][read_num])

    loggers[hpgl_id]['SL'][read_num].info("# Finished removing false hits.")

    # Let Ruffus know we are done
    open(output_file, 'w').close()

@transform(remove_false_hits,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).remove_false_hits'),
           r'\1/\2_\3.map_sl_reads',
           r'\2', r'\3')
def map_sl_reads(input_file, output_file, hpgl_id, read_num):
    """Maps the filtered spliced-leader containing reads back to the genome"""
    output_dir = '%s/%s/tophat/%s_sl_reads' % (sl_build_dir, hpgl_id, read_num)
    genome = os.path.splitext(args.genome)[0]

    loggers[hpgl_id]['SL'][read_num].info("# Mapping filtered reads back to genome")

    # input read base directory
    basedir = '%s/%s/fastq' % (sl_build_dir, hpgl_id)

    # R1 filepath (including matched SL sequence)
    if (read_num == 'R1'):
        r1_filepath = (
            '%s/actual_sl_reads/%s_R1_1_without_sl.fastq.gz' %
            (basedir, hpgl_id)
        )

        # R2 filepath (for PE reads)
        r2_filepath = r1_filepath.replace('1_with_sl', '2')

        # If SE, set filepath to empty string
        if not os.path.exists(r2_filepath):
            r2_filepath = ""
    # R2
    else:
        r2_filepath = (
            '%s/actual_sl_reads/%s_R2_2_without_sl.fastq.gz' %
            (basedir, hpgl_id)
        )
        r1_filepath = r2_filepath.replace('2_with_sl', '1')

    # Map reads using Tophat
    #  --no-mixed ?
    ret = run_tophat(output_dir, genome, loggers[hpgl_id]['SL'][read_num],
                 r1_filepath, r2_filepath,
                 extra_args='--mate-inner-dist 170 --transcriptome-max-hits 1')

    # Make sure tophat succeeded
    if ret != 0:
        logging.error(
            "# Error running tophat 2/2! %s (%s)" % (hpgl_id, read_num)
        )
        sys.exit()

    loggers[hpgl_id]['SL'][read_num].info("# Finished mapping hits to genome")

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
    input_globstr = ('%s/*/tophat/*_sl_reads/accepted_hits_sorted.bam' %
                     sl_build_dir)

    # Output: no genes nearby
    no_nearby_genes = csv.writer(open('%s/no_matches.csv' % sl_build_dir, 'w'))
    no_nearby_genes.writerow(['read_id', 'chromosome', 'strand', 'position'])

    # Output: predicted site inside a CDS
    inside_cds = csv.writer(open('%s/inside_cds.csv' % sl_build_dir, 'w'))
    inside_cds.writerow(['read_id', 'chromosome', 'strand', 'position'])

    # Itereate over mapped reads
    for filepath in glob.glob(input_globstr):
        # get hpgl id
        hpgl_id = re.match('.*(HPGL[0-9]+).*', filepath).groups()[0]

        # file to save results for individual sample
        base_dir = filepath[:re.search('tophat', filepath).start()]
        sample_hits = os.path.join(base_dir, hpgl_id + '_sl_coordinates.csv')
        sample_csv_writer = csv.writer(open(sample_hits, 'w'))
        sample_csv_writer.writerow(['read_id', 'gene_id', 'chromosome', 
                                    'strand', 'position', 'distance'])

        # keep track of how many reads were found in the expected location
        num_good = 0
        num_no_nearby_genes = 0
        num_inside_cds = 0

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

            # first, check to make sure the acceptor site does not fall within
            # a known CDS -- if it does, save to a separate file to look
            # at later
            nearby = chromosomes[chromosome][pos - 5000:pos + 5000]
            for feature in nearby.features:
                if ((feature.location.start <= 5000) and 
                    (feature.location.end >= 5000)):
                    # predicted acceptor site is inside a CDS
                    inside_cds.writerow([read.qname, chromosome, strand, pos])
                    num_inside_cds = num_inside_cds + 1
                    continue

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
                no_nearby_genes.writerow([read.qname, chromosome, strand, pos])
                num_no_nearby_genes = num_no_nearby_genes + 1
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
        loggers[hpgl_id]['SL'][read_num].info(
            "# Found %d reads with predicted acceptor site at expected location"
            % num_good)
        loggers[hpgl_id]['SL'][read_num].info(
            "# Found %d reads with predicted acceptor site inside a known CDS"
            % num_inside_cds)
        loggers[hpgl_id]['SL'][read_num].info(
            "# Found %d reads with predicted acceptor site not proximal to any CDS"
            % num_no_nearby_genes)

    # Output filepath
    coordinates_output = '%s/sl_coordinates.csv' % (sl_build_dir)
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

# run pipeline
if __name__ == "__main__":
    pipeline_run([compute_coordinates], logger=logging.getLogger(''),
                 multiprocess=args.num_threads)
    pipeline_printout_graph("utr_analysis_flowchart.png", "png",
                            [compute_coordinates])

