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

"""
import os
import re
import csv
import sys
import glob
import gzip
import time
import pysam
import random
import logging
import argparse
import datetime
import distance
import StringIO
import textwrap
import warnings
import subprocess
from ruffus import *
from Bio import Seq,SeqIO,BiopythonDeprecationWarning
from BCBio import GFF

# Hide Biopython deprecation warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

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
    ./utr_analysis.py                                               \\
        -i "$RAW/tcruzir21/*/processed/*.filtered.fastq.gz"            \\
        -s AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG                  \\
        -f1 TriTrypDB-7.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta \\
        -f2 mm10.fasta                                              \\
        -g TrypDB-7.0_TcruziCLBrenerEsmeraldo-like.gff              \\
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
                        help=('RNA-Seq FASTQ or gzipped FASTQ glob string or'
                              'a txt file containing filepaths to the samples'
                              'to be used'))
    parser.add_argument('-d', '--build-directory', required=True,
                        help='Directory to save output to')
    parser.add_argument('-f1', '--target-genome', dest='target_genome',
                        required=True, help=('Genome sequence FASTA filepath '
                        'for target species'))
    parser.add_argument('-f2', '--nontarget-genome', dest='nontarget_genome',
                        help=('Genome sequence FASTA filepath for species to '
                              'be filtered out prior to mapping. (optional)'))
    parser.add_argument('-g1', '--target-annotations', dest='target_gff',
                        required=True, help='Genome annotation GFF')
    parser.add_argument('-g2', '--nontarget-annotations', dest='nontarget_gff',
                        help='Genome annotation GFF')
    parser.add_argument('-u', '--gff-uorf-annotations', dest='uorf_gff',
                        help=('GFF containing possible uORFs; generated during '
                              'final steps of processing and can be passed '
                              'back in to re-perform analyses treating these '
                              'locations as putative ORFs.'))
    parser.add_argument('-s', '--sl-sequence', dest='spliced_leader',
                        required=True, help='Spliced leader DNA sequence',
                        default=None)
    parser.add_argument('--exclude-internal-sl-matches',
                        help=('Only allow matches with the SL at the upstream'
                              'end of a read.'), action='store_true')
    parser.add_argument('--exclude-internal-polya-matches',
                        help=('Only allow matches with the Poly(A) tail at the'
                              'downstream end of a read.'), action='store_true')
    parser.add_argument('--max-dist-from-edge',
                        help=('For unanchored searches, what is the maximum '
                        'distance from the edge of the read for a feature '
                        'match to be considered.'), default=30)
    parser.add_argument('-m', '--min-sl-length', default=10, type=int,
                        help='Minimum length of SL match (default=10)')
    parser.add_argument('-p', '--min-polya-length', default=10, type=int,
                        help='Minimum length of Poly-A match (default=10)')
    parser.add_argument('-w', '--window-size', default=15000, type=int,
                        help=('Number of bases up or downstream of feature to'
                              'scan for related genes (default=15000)'))
    parser.add_argument('-x', '--minimum-differences', default=2, type=int,
                        help=('Minimum number of differences from genomic '
                              ' sequence for a hit to be considered real. '
                              '(default=2)'))
    parser.add_argument('--num-threads', default=4, type=int,
                        help='Number of threads to use (default=4).')

    # Parse arguments
    args = parser.parse_args()

    # Replace any environmental variables and return args
    args.target_genome = os.path.expandvars(args.target_genome)
    args.input_reads = os.path.expandvars(args.input_reads)
    args.target_gff = os.path.expandvars(args.target_gff)

    if args.nontarget_genome:
        args.nontarget_genome = os.path.expandvars(args.nontarget_genome)
    if args.nontarget_gff:
        args.nontarget_gff = os.path.expandvars(args.nontarget_gff)

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
                # 2013/03/31 Bug fix (keith)
                if l[0] in '>@': # fasta/q header line
                #if l[:4] == '@HWI': # Illumina fastq header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        # Modified to preserve full identifier
        # Keith 2014-01-31
        # name, seqs, last = last[1:].partition(" ")[0], [], None
        name, seqs, last = last, [], None
        for l in fp: # read the sequence
            # Keith 2014/03/31
            if l[0] in '@+>':
            #if l[0] == '+' or l[:4] == '@HWI':
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

def run_command(cmd, log_handle, wait=True):
    """Runs a command and logs the output to a specified log handle"""
    log_handle.info("# " + cmd)

    # asynchronous
    if wait is False:
        process = subprocess.Popen(cmd.split(" "))
        return 0

    # synchronous
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

    # delete unsorted verion
    os.remove(base_output + ".bam")

    # index bam
    index_cmd = 'samtools index %s' % (base_output + '_sorted.bam')

    # 2014/04/01
    # samtools index sometimes gets stuck during execution, even after the
    # indexing has finished; since it isn't necessary for downstream processes
    # the indexing will be done asynchronously for now.
    run_command(index_cmd, log_handle, wait=False)
    log_handle.info("# Done sorting and indexing")

def run_tophat(output_dir, genome, log_handle, r1, r2="", gff=None,
               num_threads=1, read_mismatches=2, max_multihits=20, extra_args=""):
    """
    Uses Tophat to map reads with the specified settings.

    Several useful Tophat parameters have been included as keyword arguments,
    along with their default settings for Tophat 2.

    The documentation for these and other useful parameters for this script are
    included below for convenience.

    -N/--read-mismatches
        Final read alignments having more than these many mismatches are
        discarded. The default is 2.
        
    --no-novel-juncs
        Only look for reads across junctions indicated in the supplied GFF or 
        junctions file. (ignored without -G/-j)
        
    -G/--GTF <GTF/GFF3 file>	
        Supply TopHat with a set of gene model annotations and/or known
        transcripts, as a GTF 2.2 or GFF3 formatted file. If this option is
        provided, TopHat will first extract the transcript sequences and use
        Bowtie to align reads to this virtual transcriptome first. Only the
        reads that do not fully map to the transcriptome will then be mapped on
        the genome. 

    -x/--transcriptome-max-hits
        Maximum number of mappings allowed for a read, when aligned to the
        transcriptome (any reads found with more then this number of mappings
        will be discarded).

    -p/--num-threads <int>
        Use this many threads to align reads. The default is 1.

    -g/--max-multihits <int>
        Instructs TopHat to allow up to this many alignments to the reference
        for a given read, and choose the alignments based on their alignment
        scores if there are more than this number. The default is 20 for read
        mapping. Unless you use --report-secondary-alignments, TopHat will
        report the alignments with the best alignment score. If there are more
        alignments with the same score than this number, TopHat will randomly
        report only this many alignments. In case of using
        --report-secondary-alignments, TopHat will try to report alignments up
        to this option value, and TopHat may randomly output some of the
        alignments with the same score to meet this number.
    """
    # create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o755)

    # FASTA filename without extension
    genome_basename = os.path.splitext(genome)[0]

    # @NOTE 2014/02/06 -- Tophat 2.0.10 fails for some reads when attempting
    # to use more than one thread. For now, just perform the mapping
    # using a single-thread to be safe...

    # build command
    if gff is not None:
        cmd = 'tophat --no-novel-juncs -G %s ' % gff
    else:
        cmd = 'tophat '

    cmd += "--num-threads %d --max-multihits %d --read-mismatches %d %s -o %s %s %s %s" % (
           num_threads, max_multihits, read_mismatches, extra_args,
           output_dir, genome_basename, r1, r2)

    # run tophat
    ret = run_command(cmd, log_handle)

    # check to see if tophat succeeded and stop execution otherwise
    if ret != 0:
        log_handle.error("# Error running tophat (%s)!" % genome)
        print("# Error running tophat (%s)!" % genome)
        sys.exit()

    # sort and index bam output using samtools
    log_handle.info("# Sorting and indexing Tophat output")
    sort_and_index(os.path.join(output_dir, 'accepted_hits'), log_handle)
    sort_and_index(os.path.join(output_dir, 'unmapped'), log_handle)

    return 0

def run_bam2fastx(bam, fastq, log_handle):
    """
    Uses the tophat bam2fastx tool to convert a bam file to fastq

    Parameters
    ----------
    bam: str
        Input bam file
    fastq: str
        Output fastq file
    """
    #bam2fastx parameters:
    #    -q fastq
    #    -A all reads
    #    -P pair-end
    #    -o output filepath
    cmd = "bam2fastx -q -A -P -o %s %s" % (fastq, bam)

    return run_command(cmd, log_handle)

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

def filter_mapped_reads(r1, r2, genome, tophat_dir, output_fastq, log_handle,
                        gff=None, read_mismatches=2):
    """
    Maps reads using tophat and discards any that align to the genome.

    Parameters
    ----------
    r1: str
        Filepath to R1 reads to be filtered.
    r2: str
        Filepath to R2 reads to be filtered.
    genome: str
        Filepath to genome FASTA to map reads against
    tophat_dir: str
        Directory to save tophat output to
    output_fastq: str
        Filepath to save filtered FASTQ files to
    log_handle: logging.Handle
        Handler to use for logging.
    gff: string
        Filepath to GFF to use during mapping (optional)
    read_mismatches: int
        Number of mismatches to allow when scanning with Tophat.
    """
    # bam / fastq filepaths
    bam_input = os.path.join(tophat_dir, 'unmapped_sorted.bam')

    # check to see if tophat output already exists
    if not os.path.exists(bam_input):
        # map reads to genome using tophat
        log_handle.info("# Mapping against %s" % os.path.basename(genome))
        ret = run_tophat(tophat_dir, genome, log_handle, r1, r2,
                         read_mismatches=read_mismatches, gff=gff,
                         extra_args='--no-mixed')

        # number of reads before filtering
        num_reads_total = num_lines(r1) / 4
        num_unmapped = num_lines(bam_input) / 2

        # delete uneeded accepted_hits files
        os.remove(os.path.join(tophat_dir, "accepted_hits_sorted.bam"))

        log_handle.info(
            "# Ignoring %d reads which mapped to specified genome (%d remaining)" %
            (num_reads_total - num_unmapped, num_unmapped)
        )
    else:
        log_handle.info("# Skipping %s: output already exists" % tophat_dir)

    # keep unampped reads
    log_handle.info("# Converting Tophat bam output to FASTQ")
    run_bam2fastx(bam_input, output_fastq, log_handle)

def get_next_file_name(base_name):
    """Returns a filepath for the next highest file number, e.g.
       file.5.log."""
    if not os.path.exists(base_name):
        return base_name
    else:
        file_nums = [int(x.split('.').pop())
                        for x in glob.glob("%s.*" % base_name)]
        next_file_num = max([0] + file_nums) + 1
        return "%s.%d" % (base_name, next_file_num)

def find_sequence(input_file, feature_name, sequence_filter, feature_regex,
                  build_dir, sample_id, read_num, trimmed_side='left',
                  reverse=False):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing a specified sequence of interest.

    Parameters
    ----------
    input_file: str
        Filepath to a FASTQ file containing reads to scan.
    feature_name: str
        Type of feature being searched for; used in naming filing and
        directories and in choosing logs to write to. [sl|rsl|polya|polyt]
    sequence_filter: str
        A short sequence string used for initial filtering. All reads will be
        checked to see if it contains this string, and those that do will be
        further checked using a regular expression to find the location of the
        match.
    feature_regex: str
        A regular expression string indicating the exact sequence to be
        searched for. This will be either a set of spliced leader prefixes or
        suffixes, or a string of A's or T's, possibly anchored at one end of
        the read.
    build_dir: str
        Base directory to save output to.
    sample_id: str
        ID of the sample being scanned.
    read_num: str
        Which of the mated reads should be scanned. [R1|R2]
    trimmed_side: str
        When trimming the matched sequence, also remove all bases in this
        direction of the match (default: left)
    reverse: bool
        If set to True, the reverse complement of each of the trimmed reads
        will be saved.

    Output files
    ------------
    There are three possible sets of output files for this function depending
    on whether the input read comes from a mated-pair of reads or a single
    read, and whether (for the case of mated-pair reads) it is the left read
    or right read:

    1. Sequence found in R1
        *_R1_1_xxx_untrimmed.fastq
        *_R1_1_xxx_trimmed.fastq
        *_R1_2.fastq
    2. Sequence found in R2
        *_R2_2_xxx_untrimmed.fastq
        *_R2_2_xxx_trimmed.fastq
        *_R2_1.fastq
    """
    # logger
    log = loggers[sample_id][feature_name][read_num]
    log.info("# Processing %s" % os.path.basename(input_file))

    # determine whether regular expression in anchored
    anchored = (feature_regex.startswith("^") or feature_regex.endswith("$"))

    # list to keep track of potential matches
    matches = []

    # output filepaths
    output_base = '%s/%s/fastq/unfiltered/%s_%s_%s' % (
        build_dir, sample_id, sample_id, read_num, read_num[-1]
    )
    output_untrimmed = "%s_%s_untrimmed.fastq.gz" % (output_base, feature_name)
    output_trimmed = "%s_%s_trimmed.fastq.gz" % (output_base, feature_name)

    # Also keep track of match lengths which will be used for more rigorous
    # filtering when comparing to the genome sequence near where the read is
    # mapped.
    match_lengths_dir = os.path.join(build_dir, sample_id, 'results')
    output_lengths = "%s/match_lengths_%s.csv.gz" % (match_lengths_dir, read_num)
    match_lengths_fp = gzip.open(output_lengths, 'wb')

    # mated reads
    read_num_other = "R1" if read_num == "R2" else "R2"
    input_file_mated = input_file.replace("." + read_num[-1],
                                          "." + read_num_other[-1])
    output_mated_reads = "%s_%s.fastq.gz" % (output_base[:-2], read_num_other)

    # compile regex
    read_regex = re.compile(feature_regex)

    # total number of reads
    num_reads = num_lines(input_file) / 4

    # Start sample log
    log.info("# Scanning %d reads for %s" % (num_reads, feature_name))
    log.info("# Using Regex patten:\n %s" % feature_regex)

    # open output string buffer (will write to compressed file later)
    reads_trimmed = StringIO.StringIO()
    reads_untrimmed = StringIO.StringIO()
    mated_reads_buffer = StringIO.StringIO()

    # Keep track of matched read IDs
    read_ids = []

    # Find all reads containing the sequence of interest
    fastq = gzip.open(input_file, 'rb')
    fastq_mated = gzip.open(input_file_mated, 'rb')

    # iterate over mated reads at same time
    mated_reads = readfq(fastq_mated)

    for i, read in enumerate(readfq(fastq)):
        # get mated read
        mated_read = mated_reads.next()

        # ignore any reads that don't contain at least the smallest part of
        # the sequence of interest
        # this just speeds up the search so we don't have to use regex on all
        # reads
        if sequence_filter not in read[SEQUENCE_IDX]:
            continue

        # check for match

        # When looking for internal matches, there may be multiple hits. Choose
        # the one that is closest to the edge of the read where the feature is
        # expected to be found.
        try:
            # reverse sequence to search from right to left
            # in case of the reverse SL, regex has already been inverted as
            # well
            if (trimmed_side == 'right' and (not anchored)):
                match = re.search(read_regex, read[SEQUENCE_IDX][::-1])
                match_start = len(read[SEQUENCE_IDX]) - match.end() 
                match_end = len(read[SEQUENCE_IDX]) - match.start() 
            else:
                match = re.search(read_regex, read[SEQUENCE_IDX])

                # move on the next read if no match is found
                if match is None:
                    continue

                match_start = match.start()
                match_end = match.end()
        except:
            import pdb; pdb.set_trace()

        # match length
        match_length = match.end() - match.start()

        # If feature is expected to be found at left of read, trim everything
        # up to the end of the match
        if trimmed_side == 'left':
            trimmed_read = [read[ID_IDX],
                            read[SEQUENCE_IDX][match_end:],
                            "+",
                            read[QUALITY_IDX][match_end:]]
        else:
            # otherwise trim from the start of the match to the end of the read
            trimmed_read = [read[ID_IDX],
                            read[SEQUENCE_IDX][:match_start],
                            "+",
                            read[QUALITY_IDX][:match_start]]

        # skip reads that are less than 36 bases after trimming
        if len(trimmed_read[SEQUENCE_IDX]) < 36:
            continue

        # length of portion trimmed off
        trimmed_part_length = (len(read[SEQUENCE_IDX]) - 
                               len(trimmed_read[SEQUENCE_IDX]))

        # for internal matches, skip reads where match is not close enough to
        # the edge of the read
        if (trimmed_part_length - match_length) > args.max_dist_from_edge:
            continue

        # write length
        match_lengths_fp.write(",".join([read[ID_IDX], str(match_length)]) + "\n")

        # take reverse complement if requested
        # this will return the read back to the expected orientation (SL
        # upstream/Poly(A) downstream)
        #if reverse:
        #    trimmed_read[SEQUENCE_IDX] = str(
        #        Seq.Seq(trimmed_read[SEQUENCE_IDX]).reverse_complement())
        #    trimmed_read[QUALITY_IDX] = trimmed_read[QUALITY_IDX][::-1]

        # Otherwise add trimmed read to output
        reads_trimmed.write("\n".join(trimmed_read) + "\n")

        # Also save complete (untrimmed) reads containing the matched sequence.
        # By mapping these reads to the genome we can eliminate false hits;
        # i.e. reads that contain a portion of the sequence of intereste but
        # are not actual trans-splicing / poly-adenylation reads.
        untrimmed_read = [read[ID_IDX],
                          read[SEQUENCE_IDX],
                          "+",
                          read[QUALITY_IDX]]
        reads_untrimmed.write("\n".join(untrimmed_read) + "\n")

        # paired-end reads
        untrimmed_mated_read = [mated_read[ID_IDX],
                                mated_read[SEQUENCE_IDX],
                                "+",
                                mated_read[QUALITY_IDX]]
        mated_reads_buffer.write("\n".join(untrimmed_mated_read) + "\n")

        # save id
        read_ids.append(read[ID_IDX])

    # log numbers
    log.info("# Found %d reads with possible %s fragment" %
             (len(read_ids), feature_name))

    # Create output directory
    output_dir = os.path.dirname(output_base)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o755)

    # write trimmed and untrimmed reads to fastq.gz
    gzip_str(output_trimmed, reads_trimmed)
    gzip_str(output_untrimmed, reads_untrimmed)
    gzip_str(output_mated_reads, mated_reads_buffer)

    # clean up
    fastq.close()
    fastq_mated.close()
    reads_trimmed.close()
    reads_untrimmed.close()
    mated_reads_buffer.close()
    match_lengths_fp.close()

    log.info("# Finished processing %s" % os.path.basename(input_file))

def map_reads(feature_name, build_dir, sample_id, read_num):
    """
    Maps the filtered reads back to the genome. If the trimmed version of
    the read successfully maps this is indicative of a possible trans-spliced
    or polyadenylated read.
    """
    output_dir = '%s/%s/tophat/%s_filtered_trimmed' % (
        build_dir, sample_id, read_num
    )

    loggers[sample_id][feature_name][read_num].info(
        "# Mapping filtered reads back to genome"
    )

    # input read base directory
    basedir = '%s/%s/fastq' % (build_dir, sample_id)

    # R1 input filepath (including matched sequence)
    if (read_num == 'R1'):
        r1_filepath = (
            '%s/unfiltered/%s_R1_1_%s_trimmed.fastq.gz' %
            (basedir, sample_id, feature_name)
        )

        # R2 filepath (for PE reads)
        r2_filepath = r1_filepath.replace('1_%s_untrimmed' % feature_name, '2')

        # If SE, set filepath to empty string
        if not os.path.exists(r2_filepath):
            r2_filepath = ""
    # R2 input filepath
    else:
        r2_filepath = (
            '%s/unfiltered/%s_R2_2_%s_trimmed.fastq.gz' %
            (basedir, sample_id, feature_name)
        )
        r1_filepath = r2_filepath.replace('2_%s_untrimmed' % feature_name, '1')

    # Map reads using Tophat
    #  --no-mixed ?
    ret = run_tophat(output_dir, args.target_genome,
                 loggers[sample_id][feature_name][read_num],
                 r1_filepath, r2_filepath, gff=args.target_gff,
                 extra_args='--transcriptome-max-hits 1')

    loggers[sample_id][feature_name][read_num].info(
        "# Finished mapping hits to genome"
    )


def compute_coordinates(feature_name, build_dir, sample_id, read_num):
    """
    Outputs coordinates and frequencies of putative UTR features to a GFF file.

    Sequences manipulated below include:

        untrimmed_seq       untrimmed read sequence
        trimmed_seq         sequence of read after trimming
        trimmed_part     sequence of the portion of the read trimmed off
        trimmed_genome_seq  sequence at the genome location corresponding to
                            trimmed portion of the read
        matched_seq         sequence matched in previous step
        matched_genome_seq  sequence at the genome location corresponding to
                            the matched portion of the read

    Note that, unless there is unanchored matching is used and the match is
    internal (e.g. there is unrelated sequence upstream of the match), the
    trimmed_part and matched_seq will be the same, otherwise the matched
    sequences will be subsets of the trimmed ones.


    References
    ----------
    * http://www.cgat.org/~andreas/documentation/pysam/api.html#pysam.Samfile
    * http://biopython.org/wiki/GFF_Parsing
    * http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html
    """
    # output directory
    output_dir = '%s/%s/results' % (build_dir, sample_id)

    # Load genome sequence
    genome = SeqIO.parse(args.target_genome, 'fasta')
    chr_sequences = {x.id:x for x in genome}

    # Get chromosomes from GFF file
    chromosomes = load_annotations()

    # Create a dictionary to keep track of the coordinates.
    results = {}

    # keep track of how many reads were found in the expected location
    num_good = 0
    num_no_nearby_genes = 0
    num_inside_cds = 0

    # Create output CSV writer for hits that are not near any genes
    no_nearby_genes = csv.writer(
        open('%s/no_nearby_genes_%s.csv' % (output_dir, read_num), 'w')
    )
    no_nearby_genes.writerow(['read_id', 'chromosome', 'strand', 'position'])

    # Create output CSV writer for hits that are found inside known CDS's
    inside_cds = csv.writer(
        open('%s/inside_cds_%s.csv' % (output_dir, read_num), 'w')
    )
    inside_cds.writerow(['read_id', 'ch-romosome', 'strand', 'position'])

    # In addition to saving the coordinates as a GFF file, we will also write a
    # CSV file which contains entries for each read used. This can be useful
    # for debugging/tracking down the origin of a particular coordinate.
    sample_csv_writer = csv.writer(
        open('%s/matched_reads_%s.csv' % (output_dir, read_num), 'w')
    )
    sample_csv_writer.writerow(['read_id', 'gene_id', 'chromosome',
                                'strand', 'trimmed_start', 'trimmed_stop',
                                'untrimmed_start', 'untrimmed_stop',
                                'untrimmed_seq', 'trimmed_part',
                                'trimmed_genome_seq',
                                'acceptor_site'])

    # Load the untrimmed reads in order to to determine original read lengths
    input_bam_untrimmed = '%s/%s/tophat/mapped_to_target_untrimmed/unmapped_sorted.bam' % (
        shared_build_dir, sample_id
    )
    sam_untrimmed = pysam.Samfile(input_bam_untrimmed, 'rb')

    untrimmed_reads = {}

    for read in sam_untrimmed:
        rnum = 'R1' if read.is_read1 else 'R2'
        if not read.qname in untrimmed_reads:
            untrimmed_reads[read.qname] = {rnum:read}
        else:
            untrimmed_reads[read.qname][rnum] = read

    # open trimmed reads sam file
    input_bam = '%s/%s/tophat/%s_filtered_trimmed/accepted_hits_sorted.bam' % (
        build_dir, sample_id, read_num
    )
    sam = pysam.Samfile(input_bam, 'rb')

    # open match lengths file
    match_lengths_fp = gzip.open("%s/match_lengths_%s.csv.gz" % (
        output_dir, read_num), 'rb')
    match_reader = csv.reader(match_lengths_fp)

    match_lengths = {}

    for row in match_reader:
        # drop the @ in front of HWI
        match_lengths[row[0][1:]] = int(row[1])

    # Get coordinate and strand for each read in bam file
    for read in sam:
        # Get read where feature sequence was found
        if ((read.is_read1 and read_num == 'R2') or
            (read.is_read2 and read_num == 'R1')):
            continue

        # Chromosome and strand where the read was mapped
        chromosome = sam.getrname(read.tid)
        strand = "-" if read.is_reverse else "+"

        # Length of matched SL suffix or number of A's/T's
        untrimmed_read = untrimmed_reads[read.qname][read_num]

        # NOTE: this may be larger than the SL sequence length if unanchored
        # matches are allowed
        trimmed_part_length = untrimmed_read.rlen - read.rlen

        # get match length (actual portion matched)
        try:
            match_length = match_lengths[read.qname]
        except:
            import pdb; pdb.set_trace()

        # Side of original read trimmed
        if feature_name in ['sl', 'polyt']:
            trimmed_side = 'left'
        elif feature_name in ['polya', 'rsl']:
            trimmed_side = 'right'

        # Actual strand of sequence
        if feature_name in ['sl', 'polya']:
            actual_strand = strand
        elif feature_name in ['rsl', 'polyt']:
            # For rsl/polyt, the actual strand is the opposite of where read
            # was mapped
            actual_strand = '+' if strand is '-' else '-' 

        # Side of acceptor site after mapping
        if strand == '+':
            acceptor_site_side = trimmed_side
        else:
            acceptor_site_side = 'left' if trimmed_side == 'right' else 'right'

        # Acceptor site position
        if acceptor_site_side == 'left':
            acceptor_site = read.pos
        else:
            acceptor_site = read.pos + read.rlen + 1

        # Minimum length of feature
        if feature_name in ['sl', 'rsl']:
            min_feature_length = args.min_sl_length
        else:
            min_feature_length = args.min_polya_length

        # If read is mapped too close to end of chromosome, skip it
        # left end
        if read.pos < (min_feature_length - 1):
            continue
        # right end
        if ((read.pos + read.rlen) > (len(chromosomes[chromosome]) - 
             min_feature_length)):
            continue
        
        # Genome sequence just upstream of mapped location
        if acceptor_site_side == 'left': 
            # Shorten read if left end of read was mapped near end of
            # chromosome
            trimmed_part_length = min(min(trimmed_part_length, read.pos),
                                      len(chr_sequences[chromosome]) - read.pos)

            # Grab region just before mapped read
            trimmed_genome_seq = chr_sequences[chromosome][read.pos -
                trimmed_part_length:read.pos].seq

            matched_genome_seq = trimmed_genome_seq[-match_length:]

        # Genome sequence just downstream of mapped location
        else:
            # Shorten if right end of read was mapped near end of chromosome
            trimmed_part_length = min(min(trimmed_part_length, read.pos),
                                      len(chr_sequences[chromosome]) - 
                                      (read.pos + read.rlen))

            # Grab region just after mapped untrimmed read
            trimmed_genome_seq = chr_sequences[chromosome][read.pos +
                read.rlen:read.pos + read.rlen + trimmed_part_length].seq
            matched_genome_seq = trimmed_genome_seq[:match_length]

        # For reads mapped to negative strand, take complement
        if strand == '+':
            trimmed_genome_seq = str(trimmed_genome_seq)
            matched_genome_seq = str(matched_genome_seq)
        else:
            trimmed_genome_seq = str(trimmed_genome_seq.reverse_complement())
            matched_genome_seq = str(matched_genome_seq.reverse_complement())

        # Get sequence of the trimmed portion of the read and the sequence that
        # matched feature of interest
        if trimmed_side == 'left':
            # get trimmed portion from left of read
            trimmed_part = str(untrimmed_read.seq[:trimmed_part_length])
            matched_seq = trimmed_part[-match_length:]
        else:
            trimmed_part = str(untrimmed_read.seq[-trimmed_part_length:])
            matched_seq = trimmed_part[:match_length]

        # For Poly(A) tail, extend detected read to include any A's present in
        # the genome that were trimmed off
        if feature_name == 'polya':
            # + strand
            if strand == '+':
                match = re.search('^A*', matched_genome_seq)
                overlap_length = match.end() - match.start()

                if overlap_length > 0:
                    # give back some A's and update the relevant sequences
                    matched_seq = matched_seq[overlap_length:]
                    matched_genome_seq = matched_genome_seq[overlap_length:]
                    #matched_genome_seq = (matched_genome_seq + 
                    #                      ("A" * overlap_length))
                    acceptor_site = acceptor_site + overlap_length
            # - strand
            elif strand == '-':
                match = re.search('T*$', matched_genome_seq)
                overlap_length = match.end() - match.start()

                if overlap_length > 0:
                    # TESTING
                    #import pdb; pdb.set_trace()

                    # give back some A's and update the relevant sequences
                    matched_seq = matched_seq[:-overlap_length]
                    matched_genome_seq = matched_genome_seq[:-overlap_length]
                    #matched_genome_seq = (("T" * overlap_length) +
                    #                       matched_genome_seq)
                    acceptor_site = acceptor_site - overlap_length

        # Poly(T)
        elif feature_name == 'polyt':
            # + strand (corresponds to gene on negative strand)
            if strand == '+':
                match = re.search('T*$', matched_genome_seq)
                overlap_length = match.end() - match.start()

                if overlap_length > 0:
                    # TESTING
                    #import pdb; pdb.set_trace()

                    # give back some A's and update the relevant sequences
                    matched_seq = matched_seq[:-overlap_length]
                    matched_genome_seq = matched_genome_seq[:-overlap_length]
                    #matched_genome_seq = (matched_genome_seq + 
                    #                      ("T" * overlap_length))
                    acceptor_site = acceptor_site - overlap_length
            # - strand
            elif strand == '-':
                match = re.search('^A*', matched_genome_seq)
                overlap_length = match.end() - match.start()

                if overlap_length > 0:
                    # TESTING
                    #import pdb; pdb.set_trace()

                    # give back some A's and update the relevant sequences
                    matched_seq = matched_seq[overlap_length:]
                    matched_genome_seq = matched_genome_seq[overlap_length:]
                    #matched_genome_seq = (matched_genome_seq + 
                    #                      ("A" * overlap_length))
                    acceptor_site = acceptor_site + overlap_length

        # 
        # More filtering
        #

        # Check to see that match differs from genome sequence (quick)
        if (matched_seq == matched_genome_seq):
            continue

        try:
            # Make sure that for the region matched, there are at least n
            # differences between the match and the corresponding genomic
            # sequence (slower)
            match_dist = distance.hamming(matched_seq, matched_genome_seq)
            if match_dist < args.minimum_differences:
                continue
        except:
            # occurs when sl is smaller than genome/trimed_portion
            print("Trimmed: " + trimmed_part)
            print("Trimmed (genome): " + trimmed_genome_seq)
            print("Matched: " + matched_seq)
            print("Matched (genome): " + matched_genome_seq)
            import pdb; pdb.set_trace()
            continue

        # Check to make sure the acceptor site does not fall within
        # a known CDS: if it does, save to a separate file to look at later
        if is_inside_cds(chromosomes[chromosome], acceptor_site):
            inside_cds.writerow([read.qname, chromosome, strand, acceptor_site])
            num_inside_cds = num_inside_cds + 1
            continue

        # Find nearest gene
        gene = find_closest_gene(chromosomes[chromosome], feature_name,
                                 acceptor_site, acceptor_site_side)

        # If no nearby genes were found, stop here
        if gene is None:
            no_nearby_genes.writerow(
                [read.qname, chromosome, strand, acceptor_site])
            num_no_nearby_genes = num_no_nearby_genes + 1
            continue

        num_good = num_good + 1

        # Add to output dictionary
        if not chromosome in results:
            results[chromosome] = {}
        if not gene['id'] in results[chromosome]:
            results[chromosome][gene['id']] = {}

        # Add entry to sample output csv
        sample_csv_writer.writerow([
            read.qname, gene['id'], chromosome, strand,
            read.pos, read.pos + read.rlen,
            untrimmed_read.pos, untrimmed_read.pos + untrimmed_read.rlen,
            untrimmed_read.seq, trimmed_part, trimmed_genome_seq,
            acceptor_site
        ])

        # Increment site count and save distance from gene
        if not acceptor_site in results[chromosome][gene['id']]:
            results[chromosome][gene['id']][acceptor_site] = {
                "count": 1,
                "distance": gene['distance'],
                "strand": strand,
                "description": gene['description']
            }
        else:
            results[chromosome][gene['id']][acceptor_site]['count'] += 1

    # record number of good and bad reads
    loggers[sample_id][feature_name][read_num].info(
        "# Found %d reads with predicted acceptor site at expected location"
        % num_good)
    loggers[sample_id][feature_name][read_num].info(
        "# Found %d reads with predicted acceptor site inside a known CDS"
        % num_inside_cds)
    loggers[sample_id][feature_name][read_num].info(
        "# Found %d reads with predicted acceptor site not proximal to any CDS"
        % num_no_nearby_genes)

    # Output coordinates as a GFF file
    output_filepath = '%s/%s_coordinates_%s.gff' % (
        output_dir, feature_name, read_num
    )

    output_coordinates(results, feature_name, output_filepath)

    logging.info("# Finished!")

def is_inside_cds(chromosome, location):
    """
    Checks to see if a putative acceptor site is inside an annotated CDS.

    Parameters
    ----------
    chromosomes: dict
        Dictionary of chromosome SeqRecords parsed from GFF.
    location: int
        Location of putative feature of interest within the chromosome.

    Return
    ------
    bool
        True if the feature is located in a known CDS.

    References
    ----------
    http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html
    """
    # Number of bases before or after feature
    half_window = int(args.window_size / 2)

    # Extract region of sequence surrounding the putative feature of
    # interest
    nearby = chromosome[max(location - half_window, 0):location + half_window]

    # determine relative location of acceptor_site
    if location < half_window:
        # at left-end, relative site is the actual one
        relative_location = location
    else:
        # otherwise, window will be centered around the acceptor site, possibly
        # clipped at the right-end
        relative_location = len(nearby) / 2

    # scan all genes near putative feature
    for feature in nearby.features:
        if ((feature.location.start <= relative_location) and
            (feature.location.end >= relative_location)):
            return True

    return False

def find_closest_gene(chromosome, feature_name, location, acceptor_site_side):
    """
    Finds the closest gene to a specified location that is in the expected
    orientation.
    """
    # 1. Get genes within +/- args.window_size/2 bases of location
    ch_end =  int(chromosome.features[0].location.end)
    half_win = args.window_size / 2

    # Window boundaries

    # If acceptor site is on left of mapped read, look for genes to the right
    if acceptor_site_side == 'left':
        window_start = location
        window_end = min(ch_end, location + half_win)
    else:
        window_start = max(0, location - half_win)
        window_end = location

    subseq = chromosome[window_start:window_end]

    # Feature location relative to the subsequence
    feature_location = location - window_start

    # If there are no nearby genes, stop here
    if len(subseq.features) == 0:
        return None

    # Otherwise find closest gene to the acceptor site
    closest_gene = None
    closest_index = None
    closest_dist = float('inf')

    for i, gene in enumerate(subseq.features):
        # Make sure strand of gene is appropriate for the feature and
        # orientation
        if acceptor_site_side == 'left':
            if ((gene_strand == -1 and feature_name in ['sl', 'rsl']) or 
                (gene_strand ==  1 and feature_name in ['polya', 'polyt'])):
                continue
        elif acceptor_site_side == 'right':
            if ((gene_strand == -1 and feature_name in ['polya', 'polyt']) or 
                (gene_strand ==  1 and feature_name in ['sl', 'rsl'])):
                continue

        # For SL/RSL, look at gene start locations and for Poly(A)/(T) look
        # at where each gene ends.

        # Gene on positive strand
        if gene.strand == 1:
            # SL acceptor site
            if feature_name in ['sl', 'rsl']:
                dist = gene.location.start - feature_location
            # Poly(A) acceptor site
            else:
                dist = feature_location - gene.location.end

            # make sure gene is downstream of SL / upstream of Poly(A)
            if dist < 0:
                continue

        # Gene on negative strand
        else:
            # SL acceptor site
            if feature_name in ['sl', 'rsl']:
                dist = feature_location - gene.location.end
            # Poly(A) acceptor site
            else:
                dist = gene.location.start - feature_location

            # make sure gene is downstream of SL / upstream of Poly(A)
            if dist < 0:
                continue

        # For Poly(A) look at gene endings
        if dist < closest_dist:
            closest_index = i
            closest_gene = gene.id
            closest_dist = dist

    # No genes found in correct orientation
    if closest_gene is None:
        return None

    # Get description of closest matching gene
    try:
        desc_qualifiers = subseq.features[closest_index].qualifiers['description']
    except KeyError:
        # For GFFs with no description
        desc_qualifiers = []

    if len(desc_qualifiers) > 0:
        gene_description = desc_qualifiers[0]
    else:
        gene_description = ""

    # Return details for matching gene
    return {
        'id': closest_gene,
        'description': gene_description,
        'distance': closest_dist
    }

def output_coordinates(results, feature_name, filepath, track_color='0,0,255'):
    """
    Outputs the feature coordinates as a GFF3 file.

    Parameters
    ----------
    results: dict
        A nested dictionary containing the putative UTR features.
    feature_type: str
        Type of feature, as described in the SOFA ontology.
        [trans_splice_site|polyA_site]
    filepath: str
        Filepath to save results to.
    """
    fp = open(filepath, 'w')

    # Write csv header
    fp.write("##gff-version\t3\n")
    fp.write("##feature-ontology\tsofa.obo\n")
    fp.write("##attribute-ontology\tgff3_attributes.obo\n")
    fp.write("#track name=%s color=%s\n" % (feature_name.upper(), track_color))

    # Copy chromosome entries from primary GFF
    annotations_fp = open(args.target_gff)

    for line in annotations_fp:
        if "\tchromosome\t" in line:
            fp.write(line)
    annotations_fp.close()

    # Write header to output
    writer = csv.writer(fp, delimiter='\t')

    # Determine GFF feature type to use
    if feature_name in ['sl', 'rsl']: 
        feature_type = 'trans_splice_site'
    else:
        feature_type = 'polyA_site'

    # write output to csv
    for chrnum in results:
        for gene_id in results[chrnum]:
            for i, acceptor_site in enumerate(results[chrnum][gene_id], 1):
                # gff3 attributes
                attributes = "ID=%s.%s.%d;Name=%s;description=%s" % (
                    gene_id, feature_name, i, gene_id,
                    results[chrnum][gene_id][acceptor_site]['description']
                )

                # write entry
                writer.writerow([
                    chrnum, "utr_analysis.py", feature_type,
                    acceptor_site, acceptor_site,
                    results[chrnum][gene_id][acceptor_site]['count'],
                    results[chrnum][gene_id][acceptor_site]['strand'],
                    '.', attributes
                ])
    fp.close()

def output_unannotated_orfs(results, filepath):
    """
    Outputs a set of unannotated ORFs to a specified GFF file.

    Parameters
    ----------
    results: dict
        A nested dictionary containing the putative ORFs
    filepath: str
        Filepath to save results to.
    """
    fp = open(filepath, 'w')

    # Write csv header
    fp.write("##gff-version\t3\n")
    fp.write("##feature-ontology\tsofa.obo\n")
    fp.write("##attribute-ontology\tgff3_attributes.obo\n")
    fp.write("#track name=unannotated_ORFs color='134,166,83'\n")

    # Copy chromosome entries from primary GFF
    #annotations_fp = open(args.target_gff)

    #for line in annotations_fp:
    #    if "\tchromosome\t" in line:
    #        fp.write(line)
    #annotations_fp.close()

    # Create a CSV writer instance
    writer = csv.writer(fp, delimiter='\t')

    # Increment new ORFs
    i = 0

    # write output to csv
    for chrnum in results:
        for orf in results[chrnum]:
            # gff3 attributes
            orf_id = "ORF.%d" % i
            attributes = "ID=%s;Name=%s;description=%s" % (
                orf_id, orf_id, orf_id
            )

            # write entry
            writer.writerow([
                chrnum, "utr_analysis.py", "ORF",
                orf[0], orf[1], '.', orf[2], '.', attributes
            ])
            i = i + 1
    fp.close()

def combine_gff_results(input_gffs):
    """Combines GFF output generated for each individual sample into a single
    combined file."""
    # Create a dictionary to keep track of the coordinates.
    results = {}

    # GFF entry fields
    gff_fields = ['chromosome', 'script', 'feature', 'start', 'stop', 'count',
                  'strand', 'quality', 'description']

    # Parse individual GFFs and create a new summary GFF
    for gff in input_gffs:
        # Open GFF and skip first few lines (comments)
        fp = open(gff)
        fp.readline()
        fp.readline()
        fp.readline()
        fp.readline()

        # Parse with CSV reader
        reader = csv.DictReader(fp, delimiter='\t', fieldnames=gff_fields)

        for row in reader:
            # Skip chromosome entries
            if row['feature'] == 'chromosome':
                continue

            # Chromosome
            chromosome = row['chromosome']
            acceptor_site = row['start']

            count = int(row['count'])
            strand = row['strand']

            # Parse GFF attributes field
            # Example: ID=LmjF.31.ncRNA1.sl.8;Name=LmjF.31.001;description=unspecified
            parts = row['description'].split(';')

            # Parse gene name and description from attributes
            gene = parts[1][5:]
            description = parts[-1][12:]

            # Add to output dictionary
            if not chromosome in results:
                results[chromosome] = {}
            if not gene in results[chromosome]:
                results[chromosome][gene] = {}

            # Increment site count and save distance from gene
            if not acceptor_site in results[chromosome][gene]:
                results[chromosome][gene][acceptor_site] = {
                    "count": count,
                    "strand": strand,
                    "description": description
                }
            else:
                results[chromosome][gene][acceptor_site]['count'] += count

    # Return combined results
    return results

def load_annotations():
    """Loads genome annotations from specified GFF(s)."""
    # Get chromosomes from GFF file
    chromosomes = {}

    # Load existing gene annotations
    annotations_fp = open(args.target_gff)

    for entry in GFF.parse(annotations_fp):
        if len(entry.features) > 0 and entry.features[0].type == 'chromosome':
            chromosomes[entry.id] = entry

    # If file containing additional unannotated ORFs was specified, include
    # those as well
    if args.uorf_gff:
        for ch in GFF.parse(open(args.uorf_gff)):
            for feature in ch.features:
                if feature.type == 'ORF':
                    chromosomes[ch.id].features.append(feature)

    # clean up
    annotations_fp.close()

    return chromosomes

#--------------------------------------
# Main
#--------------------------------------

# parse input
args = parse_input()

# Shared directory
shared_build_dir = os.path.join(args.build_directory, 'common')

# Combined output directory
combined_output_dir = os.path.join(
    args.build_directory,
    'results',
    'sl-min%d%s_polya-min%d%s_mindiff-%d' % (
        args.min_sl_length,
        '-anchored' if args.exclude_internal_sl_matches else '',
        args.min_polya_length,
        '-anchored' if args.exclude_internal_polya_matches else '',
        args.minimum_differences
    )
)
if not os.path.exists(combined_output_dir):
    os.makedirs(combined_output_dir, mode=0o755)

# Create a unique build paths for specified parameters

# SL sub-directory
sl_build_dir = os.path.join(
    args.build_directory,
    'spliced_leader',
    'minlength%d%s_mindiff-%d' % (
        args.min_sl_length,
        '-anchored' if args.exclude_internal_sl_matches else '',
        args.minimum_differences
    )
)

# Reverse SL sub-directory
rsl_build_dir = os.path.join(
    args.build_directory,
    'reverse_spliced_leader',
    'minlength%d%s_mindiff-%d' % (
        args.min_sl_length,
        '-anchored' if args.exclude_internal_sl_matches else '',
        args.minimum_differences
    )
)

# Poly(A) tail sub-directory
polya_build_dir = os.path.join(
    args.build_directory,
    'poly-a',
    'minlength%d%s_mindiff-%d' % (
        args.min_polya_length,
        '-anchored' if args.exclude_internal_polya_matches else '',
        args.minimum_differences
    )
)

# Poly(A) tail reverse complement sub-directory
polyt_build_dir = os.path.join(
    args.build_directory,
    'poly-t',
    'minlength%d%s_mindiff-%d' % (
        args.min_polya_length,
        '-anchored' if args.exclude_internal_polya_matches else '',
        args.minimum_differences
    )
)

# Create a reversed version of the SL sequence
reverse_sl = str(Seq.Seq(args.spliced_leader).reverse_complement())

# Get a list of sample ids ids
# @TODO: Generalize handling of sample ids
input_regex = re.compile(r'.*(HPGL[0-9]+).*')
sample_ids = []

# Get samples to be parsed
filenames = glob.glob(args.input_reads)
for filename in filenames:
    sample_id = re.match(input_regex, filename).groups()[0]
    if sample_id not in sample_ids:
        sample_ids.append(sample_id)

# list of ids including previously processed samples (used for final step)
sample_ids_all = sample_ids
input_globstr = (
    '%s/*/tophat/*_sl_reads/accepted_hits_sorted.bam' % sl_build_dir
)
for filepath in glob.glob(input_globstr):
    # get hpgl id
    sample_ids_all.append(re.match('.*(HPGL[0-9]+).*', filepath).groups()[0])


# create subdirs based on matching parameters
shared_subdirs = ['ruffus',
                  'fastq/genomic_reads_removed',
                  'tophat/mapped_to_target_untrimmed']

if args.nontarget_genome:
    shared_subdirs = shared_subdirs + ['fastq/nontarget_reads_removed',
                                       'tophat/mapped_to_nontarget']

for sample_id in sample_ids:
    # shared directories
    for sub_dir in shared_subdirs:
        outdir = os.path.join(shared_build_dir, sample_id, sub_dir)
        if not os.path.exists(outdir):
            os.makedirs(outdir, mode=0o755)

    # parameter- and feature-specific directories
    for base_dir in [sl_build_dir, rsl_build_dir, 
                     polya_build_dir, polyt_build_dir]:
        for sub_dir in ['fastq/filtered', 'fastq/unfiltered',
                        'results', 'ruffus', 'log', 'tophat']:
            outdir = os.path.join(base_dir, sample_id, sub_dir)
            if not os.path.exists(outdir):
                os.makedirs(outdir, mode=0o755)

# setup master logger
log_format = '%(asctime)s %(message)s'
date_format = '%Y-%m-%d %I:%M:%S %p'
formatter = logging.Formatter(log_format, datefmt=date_format)

# determine log name to use
master_log = get_next_file_name(os.path.join(args.build_directory, 'build.log'))

logging.basicConfig(filename=master_log,
                    level=logging.INFO,
                    format=log_format,
                    datefmt=date_format)

# log to console as well
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)

# get version information
from ruffus import __version__ as ruffus_version
from Bio import __version__ as biopython_version

logging.info("# Starting UTR Analysis")
logging.info("# Python %s" % sys.version)
logging.info("# Ruffus %s" % ruffus_version)
logging.info("# Biopython %s" % biopython_version)
logging.info("# Command:\n%s" % " ".join(sys.argv))

# create dictionary of log handlers for sample-specific info
loggers = {}

# setup sample-specific loggers
for sample_id in sample_ids_all:
    loggers[sample_id] = {}

    for analysis in ['sl', 'rsl', 'polya', 'polyt']:
        loggers[sample_id][analysis] = {}

        for read_num in ['R1', 'R2']:
            if analysis == 'sl':
                bdir = sl_build_dir
            elif analysis == 'rsl':
                bdir = rsl_build_dir
            elif analysis == 'polya':
                bdir = polya_build_dir
            else:
                bdir = polyt_build_dir

            sample_log_name = get_next_file_name(
                os.path.join(bdir, sample_id, 'log', '%s_%s_%s.log' % (
                    sample_id, analysis, read_num
                ))
            )
            loggers[sample_id][analysis][read_num] = logging.getLogger(
                sample_id + analysis + read_num
            )
            handler = logging.FileHandler(sample_log_name)
            handler.setFormatter(formatter)
            loggers[sample_id][analysis][read_num].addHandler(handler)

#-----------------------------------------------------------------------------
# RUFFUS TASKS
#
# Overview
# --------
# Below are the ruffus tasks associated with the UTR analysis pipeline. The
# first few tasks (0-2) are "common" tasks that are done regardless of the
# feature of interest (SL or Poly(A)). For example -- mapping against a
# "nontarget" (e.g. host) genome to remove unrelated reads, and removing reads
# which map to the genome as-is, thus indicating they are not related to a
# splicing or polyadenylation event.
#
# Steps 3 and above are specific to the feature of interest. In each case, all
# of the steps of analysis are handled in order for each feature before moving
# on to the next feature: splice acceptor sites are process, then Poly(A) and
# then Poly(T). Each task process one read (R1 or R2) of one sample at a time.
# This is done so that multiple samples can be processed in parallel.
#
# At the very end, there is another common step which is only executed once all
# of the previous tasks have been completed. In this step, the results are
# summarized across multiple samples.
#
# Regular expressions
# -------------------
# To keep track of Ruffus's progress, empty sentinel files are created in
# the ruffus subdirectory of each sample. The sentinel filenames are encoded
# with some information about the currently running task, which are parsed
# using a regular expression. For the first main task (parse_sl_reads), the
# input is not a sentinel file, but an input fastq filepath. The various
# componenets of the regular expression used in this case is provided as
# an example below.
#
# Ex. "$RAW/tcruzir21/HPGL0121/processed/HPGL0121_R1_filtered.fastq"
#
#    \1 - directory
#    \2 - HPGLxxxx
#    \3 - _anything_between_id_and_read_num_
#    \4 - R1/R2
#    \5 - _anything_after_read_num_
#
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Step 0A: Check for bowtie indices
#
# Before beginning analysis, checks to make sure bowtie2 indices exist for the
# genomes of interest, and generates them otherwise.
#-----------------------------------------------------------------------------
def check_for_bowtie_indices():
    """check for bowtie 2 indices and create if needed"""
    # check index for target species
    genome1= os.path.splitext(args.target_genome)[0]

    # if index exists, stop here
    if os.path.exists('%s.1.bt2' % genome1):
        return

    # otherwise, create bowtie2 index
    logging.info("# Building bowtie2 index for %s" % args.target_genome)
    bowtie_cmd = (['bowtie2-build', args.target_genome, genome1])
    logging.info("# Command:\n" + " ".join(bowtie_cmd))
    ret = subprocess.call(bowtie_cmd)

    # stop here if we are only mapping to one genome
    if not args.nontarget_genome:
        return

    # check index for nontarget species
    genome2 = os.path.splitext(args.nontarget_genome)[0]

    # if index exists, stop here
    if os.path.exists('%s.1.bt2' % genome2):
        return

    # otherwise, create bowtie2 index
    logging.info("# Building bowtie2 index for %s" % args.nontarget_genome)
    bowtie_cmd = (['bowtie2-build', args.nontarget_genome, genome2])
    logging.info("# Command:\n" + " ".join(bowtie_cmd))
    ret = subprocess.call(bowtie_cmd)

#-----------------------------------------------------------------------------
# Step 0B: Check genome FASTA files
#
# Checks to make sure that the specified FASTA files exist, and that a version
# of the files ending in .fa (as required by Tophat) exist. If not, a symlinked
# version is created.
#-----------------------------------------------------------------------------
def check_genome_fastas():
    """Checks to make sure the genome fasta exists and is available as a .fa
    file; required by Tophat."""
    # First, check target genome
    target_genome = os.path.splitext(args.target_genome)[0]

    # if index exists, continue to next genome
    if os.path.exists('%s.fa' % target_genome):
        pass
    elif os.path.exists('%s.fasta' % target_genome):
        # if .fasta file exists, but not .fa, create a symlink
        logging.info("# Creating symlink to %s for Tophat" % args.target_genome)
        os.symlink(target_genome + '.fasta', target_genome + '.fa')
    else:
        raise IOError("Missing target genome file")

    # stop here if we are only mapping to one genome
    if not args.nontarget_genome:
        return

    # Next, check filter genome
    nontarget_genome = os.path.splitext(args.nontarget_genome)[0]

    # if index exists, continue to next genome
    if os.path.exists('%s.fa' % nontarget_genome):
        pass
    elif os.path.exists('%s.fasta' % nontarget_genome):
        # if .fasta file exists, but not .fa, create a symlink
        logging.info("# Creating symlink to %s for Tophat" % args.nontarget_genome)
        os.symlink(nontarget_genome + '.fasta', nontarget_genome + '.fa')
    else:
        raise IOError("Missing filter genome file")

#-----------------------------------------------------------------------------
# Step 1: Filter nontarget reads (Optional)
#
# Before looking for spliced leader (SL) and Poly(A) reads, we will first
# remove any reads which come from an unrelated species such as host cells used
# to culture the target species.
#-----------------------------------------------------------------------------
@follows(check_for_bowtie_indices)
@follows(check_genome_fastas)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq(\.gz)?'),
           r'%s/\2/ruffus/\2_\4.filter_nontarget_reads' % shared_build_dir,
           r'\2', r'\4')
def filter_nontarget_reads(input_file, output_file, sample_id, read_num):
    # If we are only mapping to a single genome, skip this step
    if not args.nontarget_genome:
        open(output_file, 'w').close()
        return

    # We only need to map once for each mated pair
    if read_num == "R2":
        # Wait for R1 task to finish processing and then mark as finished
        while not os.path.exists(output_file.replace("R2", "R1")):
            time.sleep(120)

        # Mark as finished and exit
        open(output_file, 'w').close()
        return

    logging.info("# Removing nontarget reads.")

    # output directories
    tophat_dir = os.path.join(shared_build_dir, sample_id,
                          'tophat', 'mapped_to_nontarget')
    fastq_dir = os.path.join(shared_build_dir, sample_id,
                          'fastq', 'nontarget_reads_removed')

    # read locations
    r1 = input_file
    r2 = r1.replace("R1", "R2")

    # output fastq filepaths
    output_fastq = os.path.join(
        fastq_dir, "%s_nontarget_reads_removed.fastq.gz" % (sample_id))

    # check for gff
    gff=args.nontarget_gff if args.nontarget_gff else None

    # map reads and remove hits
    filter_mapped_reads(r1, r2, args.nontarget_genome, tophat_dir,
                        output_fastq, logging, gff=gff)
    logging.info("# Finished removing nontarget reads.")

    # Let ruffus know we are finished
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Step 2: Filter genomic reads
#
# The next step is to remove reads whichS map to the target genome before any
# trimming. These reads likely represent actual features in the genome and not
# the capped or polyadenylated reads we are interested in.
#-----------------------------------------------------------------------------
@follows(filter_nontarget_reads)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq(\.gz)?'),
           r'%s/\2/ruffus/\2_\4.filter_genomic_reads' % shared_build_dir,
           r'\2', r'\4')
def filter_genomic_reads(input_file, output_file, sample_id, read_num):
    # We only need to map once for each mated pair
    if read_num == "R2":
        # Wait for R1 task to finish processing and then mark as finished
        while not os.path.exists(output_file.replace("R2", "R1")):
            time.sleep(120)

        # Mark as finished and exit
        open(output_file, 'w').close()
        return

    logging.info("# Removing genomic reads.")

    # determine source of input reads to use
    if args.nontarget_genome:
        # if nontarget genome was specified, use output from that step
        input_fastq_dir = os.path.join(shared_build_dir, sample_id,
                                    'fastq', 'nontarget_reads_removed')
        r1 = os.path.join(input_fastq_dir,
                        "%s_nontarget_reads_removed.1.fastq.gz" % (sample_id))
        r2 = r1.replace(".1", ".2")
    else:
        # otherwise use inputs specified at run-time
        r1 = input_file
        r2 = r1.replace("R1", "R2")

    # output tophat and fastq directories
    tophat_dir = os.path.join(shared_build_dir, sample_id,
                              'tophat', 'mapped_to_target_untrimmed')
    output_fastq_dir = os.path.join(shared_build_dir, sample_id,
                                   'fastq', 'genomic_reads_removed')

    # output fastq filepaths
    output_fastq = os.path.join(
        output_fastq_dir, "%s_genomic_reads_removed.fastq.gz" % (sample_id))

    # map reads and remove hits
    # Initially we will keep all (unmapped) reads which differ from genome by
    # at least one base. Later on we can be more restictive in our filtering
    # to make sure we aren't getting spurious hits.
    filter_mapped_reads(r1, r2, args.target_genome, tophat_dir, output_fastq,
                        logging, read_mismatches=1, gff=args.target_gff)
    logging.info("# Finished removing genomic reads.")

    # Let ruffus know we are finished
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Step 3: Find reads with sequence of interest
#
# Finds reads containing a minimum number of bases of the feature of interest;
# either a portion of the spliced leader (SL) sequence, or a Poly(A) or Poly(T)
# tract in the expected location for a polyadenylation event.
#
# To speed things up, the matching is done in two steps:
#
#  (1) First a simple string filter (e.g. 10 bases of the SL sequence) is used
#      to match all reads which have the minimum required similarity.
#  (2) Next, a regular expression is used to match determine the precise
#      boundaries of the match, and thus, where to trim reads.
#-----------------------------------------------------------------------------
@transform(filter_genomic_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).filter_genomic_reads'),
           r'%s/\2/ruffus/\2_\3.find_sl_reads' % sl_build_dir,
           r'\2', r'\3')
def find_sl_reads(input_file, output_file, sample_id, read_num):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing a subsequence of the spliced leader.
    """
    # First-stage filter string
    sl_filter = args.spliced_leader[-args.min_sl_length:]

    logging.info("# Finding SL reads")

    # Determine strings to match in reads
    # Note that "|" separated regexes are greedily matched from left-to-right.
    # In order to find the longest matching subsequence then we must order the
    # regex so that the longest substring comes first.
    if args.exclude_internal_sl_matches:
        subseqs = ["^" + args.spliced_leader[-x:] for x in
                  range(args.min_sl_length, len(args.spliced_leader) + 1)]
        subseqs.reverse()
        sl_regex = '|'.join(subseqs)
    else:
        subseqs = [args.spliced_leader[-x:] for x in
                   range(args.min_sl_length, len(args.spliced_leader) + 1)]
        subseqs.reverse()
        sl_regex = '|'.join(subseqs)

    # input reads
    input_reads = os.path.join(
        shared_build_dir, sample_id, 'fastq', 'genomic_reads_removed',
        "%s_genomic_reads_removed.%s.fastq.gz" % (sample_id, read_num[-1]))

    find_sequence(input_reads, 'sl', sl_filter, sl_regex, sl_build_dir,
                  sample_id, read_num)

    # Let Ruffus know we are done
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Step 4: Map trimmed reads
#
# Trims the matched sequences from reads and map to genome. For reads where
# the matched sequence comes from a trans-splicing or polyadenylation event,
# the location of the mapped trimmed read is where the addition took place.
#-----------------------------------------------------------------------------
@transform(find_sl_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).find_sl_reads'),
           r'\1/\2_\3.map_sl_reads',
           r'\2', r'\3')
def map_sl_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered spliced-leader containing reads back to the genome"""
    map_reads('sl', sl_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Step 5: Compute UTR coordinates
#
# In this step, the locations of mapped trimmed reads above are used to
# determine likely trans-splicing and polyadenylation acceptor sites in the
# genome for the processed input samples. A couple of additional filtering
# steps are based to ensure that the putative acceptor sites lie in acceptable
# locations (e.g. not inside a CDS.)
#-----------------------------------------------------------------------------
@transform(map_sl_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).map_sl_reads'),
           r'\1/\2_\3.compute_sl_coordinates',
           r'\2', r'\3')
def compute_sl_coordinates(input_file, output_file, sample_id, read_num):
    logging.info("# Computing coordinates for mapped trans-splicing events")
    compute_coordinates('sl', sl_build_dir, sample_id, read_num)
    open(output_file, 'w').close()


#-----------------------------------------------------------------------------
# Reverse SL analysis
#-----------------------------------------------------------------------------

#
# Reverse SL Step 1
#
@follows(compute_sl_coordinates)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq(\.gz)?'),
           r'%s/\2/ruffus/\2_\4.find_rsl_reads' % rsl_build_dir,
           r'\2', r'\4')
def find_rsl_reads(input_file, output_file, sample_id, read_num):
    """Matches reads with possible Poly(A) tail fragment"""
    # First-stage filter string
    rsl_filter = reverse_sl[:args.min_sl_length]

    # Determine strings to match in reads
    if args.exclude_internal_sl_matches:
        subseqs = [reverse_sl[:x] + '$' for x in
                  range(args.min_sl_length, len(args.spliced_leader) + 1)]
        subseqs.reverse()
        rsl_regex = '|'.join(subseqs)
    else:
        # In order to simplify detection of the right-most features in cases
        # where  multiple matches are present in read, we will reverse the 
        # read and regex and search and then invert results
        subseqs = [reverse_sl[:x][::-1] for x in
                   range(args.min_sl_length, len(args.spliced_leader) + 1)]
        subseqs.reverse()
        rsl_regex = '|'.join(subseqs)

    # input reads
    input_reads = os.path.join(
        shared_build_dir, sample_id, 'fastq', 'genomic_reads_removed',
        "%s_genomic_reads_removed.%s.fastq.gz" % (sample_id, read_num[-1]))

    logging.info("# find_rsl_reads (%s %s)" % (sample_id, read_num))

    find_sequence(input_reads, 'rsl', rsl_filter, rsl_regex,
                  rsl_build_dir, sample_id, read_num, trimmed_side='right',
                  reverse=True)
    open(output_file, 'w').close()

#
# RSL Step 2
#
@transform(find_rsl_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).find_rsl_reads'),
           r'\1/\2_\3.map_rsl_reads',
           r'\2', r'\3')
def map_rsl_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered poly-adenylated reads back to the genome"""
    map_reads('rsl', rsl_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#
# RSL Step 3
#
@transform(map_rsl_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).map_rsl_reads'),
           r'\1/\2_\3.compute_rsl_coordinates',
           r'\2', r'\3')
def compute_rsl_coordinates(input_file, output_file, sample_id, read_num):
    logging.info(
        "# Computing coordinates for mapped reverse sl events [A]")
    compute_coordinates('rsl', rsl_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Poly(A) Analysis
#-----------------------------------------------------------------------------

#
# Poly(A) Step 1
#
@follows(compute_rsl_coordinates)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq(\.gz)?'),
           r'%s/\2/ruffus/\2_\4.find_polya_reads' % polya_build_dir,
           r'\2', r'\4')
def find_polya_reads(input_file, output_file, sample_id, read_num):
    """Matches reads with possible Poly(A) tail fragment"""
    polya_filter = 'A' * args.min_polya_length

    if args.exclude_internal_polya_matches:
        polya_regex = 'A{%d,}$' % (args.min_polya_length)
    else:
        polya_regex = 'A{%d,}' % (args.min_polya_length)

    # input reads
    input_reads = os.path.join(
        shared_build_dir, sample_id, 'fastq', 'genomic_reads_removed',
        "%s_genomic_reads_removed.%s.fastq.gz" % (sample_id, read_num[-1]))

    logging.info("# find_polya_reads (%s %s)" % (sample_id, read_num))

    find_sequence(input_reads, 'polya', polya_filter, polya_regex,
                  polya_build_dir, sample_id, read_num, trimmed_side='right')
    open(output_file, 'w').close()

#
# Poly(A) Step 2
#
@transform(find_polya_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).find_polya_reads'),
           r'\1/\2_\3.map_polya_reads',
           r'\2', r'\3')
def map_polya_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered poly-adenylated reads back to the genome"""
    map_reads('polya', polya_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#
# Poly(A) Step 3
#
@transform(map_polya_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).map_polya_reads'),
           r'\1/\2_\3.compute_polya_coordinates',
           r'\2', r'\3')
def compute_polya_coordinates(input_file, output_file, sample_id, read_num):
    logging.info(
        "# Computing coordinates for mapped polyadenylation events [A]")
    compute_coordinates('polya', polya_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Poly(T) Analysis
#-----------------------------------------------------------------------------

#
# Poly(T) Step 1
#
@follows(compute_polya_coordinates)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq(\.gz)?'),
           r'%s/\2/ruffus/\2_\4.find_polyt_reads' % polyt_build_dir,
           r'\2', r'\4')
def find_polyt_reads(input_file, output_file, sample_id, read_num):
    """Matches reads with possible Poly(T) tail fragment"""
    # Match reads with at least n T's at the beginning of the read; For now
    # we will always require matches to be at the beginning of the read.
    polyt_filter = 'T' * args.min_polya_length
    if args.exclude_internal_polya_matches:
        polyt_regex = '^T{%d,}' % (args.min_polya_length)
    else:
        polyt_regex = 'T{%d,}' % (args.min_polya_length)

    # input reads
    input_reads = os.path.join(
        shared_build_dir, sample_id, 'fastq', 'genomic_reads_removed',
        "%s_genomic_reads_removed.%s.fastq.gz" % (sample_id, read_num[-1]))

    logging.info("# find_polyt_reads (%s %s)" % (sample_id, read_num))

    find_sequence(input_reads, 'polyt', polyt_filter, polyt_regex,
                  polyt_build_dir, sample_id, read_num, reverse=True)
    open(output_file, 'w').close()

#
# Poly(T) Step 2
#
@transform(find_polyt_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).find_polyt_reads'),
           r'\1/\2_\3.map_polyt_reads',
           r'\2', r'\3')
def map_polyt_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered Poly(T) reads back to the genome"""
    map_reads('polyt', polyt_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#
# Poly(T) Step 3
#
@transform(map_polyt_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).map_polyt_reads'),
           r'\1/\2_\3.compute_polyt_coordinates',
           r'\2', r'\3')
def compute_polyt_coordinates(input_file, output_file, sample_id, read_num):
    logging.info(
        "# Computing coordinates for mapped polyadenylation events [T]")
    compute_coordinates('polyt', polyt_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Step 6: Combine coordinate output for multiple samples
#
# Next, we create a single GFF file using the information contained in the GFF
# files that were constructed for each sample.
#
# Coordinates with low coverage may also be filtered out at this step.
#-----------------------------------------------------------------------------
@merge(compute_polyt_coordinates, '%s/finished' % combined_output_dir)
def combine_results(input_files, output_file):
    # Convert input ruffus tasks to corresponding GFF filepaths
    regex = '.*/(HPGL[0-9]+)_(R[1-2]).*'

    # Combine spliced leader output
    logging.info("# Combining spliced leader coordinate output and filtering "
                 "out sites with low support")

    sl_gffs = []
    for infile in input_files:
        (sample_id, read_num) = re.match(regex, infile).groups()
        # SL
        gff1 = "%s/%s/results/sl_coordinates_%s.gff" % (
                 sl_build_dir, sample_id, read_num)
        sl_gffs.append(gff1)

        # Reverse SL
        gff2 = "%s/%s/results/rsl_coordinates_%s.gff" % (
                 rsl_build_dir, sample_id, read_num)
        sl_gffs.append(gff2)

    sl_outfile = os.path.join(combined_output_dir, 'spliced_leader.gff')
    sl_combined = combine_gff_results(sl_gffs)

    # Combine Poly(A) output
    logging.info("# Combining Poly(A) coordinate output and filtering "
                 "out sites with low support")

    polya_gffs = []
    for infile in input_files:
        (sample_id, read_num) = re.match(regex, infile).groups()
        # Poly(A)
        gff1 = "%s/%s/results/polya_coordinates_%s.gff" % (
                 polya_build_dir, sample_id, read_num)
        polya_gffs.append(gff1)

        # Poly(T)
        gff2 = "%s/%s/results/polyt_coordinates_%s.gff" % (
                 polyt_build_dir, sample_id, read_num)
        polya_gffs.append(gff2)

    polya_outfile = os.path.join(combined_output_dir, 'polya.gff')
    polya_combined = combine_gff_results(polya_gffs)

    # Save summary GFFs
    output_coordinates(sl_combined, 'sl', sl_outfile, track_color='83,166,156')
    output_coordinates(polya_combined, 'polya', polya_outfile,
                       track_color='166,83,93')

#-----------------------------------------------------------------------------
# Run pipeline
#-----------------------------------------------------------------------------
if __name__ == "__main__":
    pipeline_run([combine_results], forcedtorun_tasks=[combine_results],
                 logger=logging.getLogger(''), multiprocess=args.num_threads)
    pipeline_printout_graph("utr_analysis_flowchart.png", "png",
                            [combine_results])

