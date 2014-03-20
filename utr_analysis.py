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
                        help=('Number of bases up or downstream of feature to'
                              'scan for related genes (default=10000)'))
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

    Several useful Tophat parameters have been included as keyword arguments,
    along with their default settings for Tophat 2.

    The documentation for these and other useful parameters for this script are
    included below for convenience.

    -N/--read-mismatches
        Final read alignments having more than these many mismatches are 
        discarded. The default is 2.

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

    # @NOTE 2014/02/06 -- Tophat 2.0.10 fails for some reads when attempting
    # to use more than one thread. For now, just perform the mapping
    # using a single-thread to be safe...

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

def find_sequence(input_file, feature_name, sequence_filter, feature_regex, 
                  build_dir, sample_id, read_num):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing a specified sequence of interest.

    Parameters
    ----------
    input_file: str
        Filepath to a FASTQ file containing reads to scan.
    feature_name: str
        Type of feature being searched for; used in naming filing and
        directories and in choosing logs to write to. [sl|polya|polyt]
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

    Output files
    ------------
    There are three possible sets of output files for this function depending
    on whether the input read comes from a mated-pair of reads or a single
    read, and whether (for the case of mated-pair reads) it is the left read
    or right read:

    1. Sequence found in R1
        *_R1_1_with_xxx.fastq
        *_R1_1_without_xxx.fastq
        *_R1_2.fastq
    2. Sequence found in R2
        *_R2_2_with_xxx.fastq
        *_R2_2_without_xxx.fastq
        *_R2_1.fastq
    """
    # sample logger
    log = loggers[sample_id][feature_name][read_num]

    # list to keep track of potential matches
    matches = []

    # output filepaths
    output_base = '%s/%s/fastq/unfiltered/%s_%s' % (
        build_dir, sample_id, sample_id, read_num 
    )
    output_with_feature = "%s_%s_with_%s.fastq" % (
        output_base, read_num[-1], feature_name
    )
    output_without_feature = "%s_%s_without_%s.fastq" % (
        output_base, read_num[-1], feature_name
    )

    # mated reads
    read_num_other = "R1" if read_num == "R2" else "R2"
    input_file_mated = input_file.replace(read_num, read_num_other)
    output_mated_reads = "%s_%s.fastq" % (output_base, read_num_other[-1])

    # compile regex
    read_regex = re.compile(feature_regex)

    # total number of reads
    num_reads = num_lines(input_file) / 4

    # Start sample log
    log.info("# Processing %s" % os.path.basename(input_file))
    log.info("# Scanning %d reads for %s" % (num_reads, feature_name))
    log.info("# Using Regex patten:\n %s" % feature_regex)

    # open output string buffer (will write to compressed file later)
    reads_without_feature = StringIO.StringIO()
    reads_with_feature = StringIO.StringIO()
    mated_reads_buffer = StringIO.StringIO()

    # Keep track of matched read IDs
    read_ids = []

    # Find all reads containing the sequence of interest
    fastq = open(input_file)
    fastq_mated = open(input_file_mated)

    # iterate over mated reads at same time
    mated_reads = readfq(fastq_mated)

    for i, read in enumerate(readfq(fastq)):
        # get mated read
        mated_read = mated_reads.next()

        # ignore any reads that don't contain at least the smallest part of
        # the sequence of interest
        if sequence_filter not in read[SEQUENCE_IDX]:
            continue

        # check for match
        match = re.search(read_regex, read[SEQUENCE_IDX])

        # move on the next read if no match is found
        if match is None:
            continue

        # otherwise add to output fastq

        # If matched sequence is at the beginning of the read, trim everything
        # up to the end of the match
        if match.start() <= (len(read[SEQUENCE_IDX]) - match.end()):
            trimmed_read = [read[ID_IDX],
                            read[SEQUENCE_IDX][match.end():],
                            "+",
                            read[QUALITY_IDX][match.end():]]
        else:
            # If matched sequence is at the end of the read, trim everything 
            # from the start of the match on
            trimmed_read = [read[ID_IDX],
                            read[SEQUENCE_IDX][:match.start()],
                            "+",
                            read[QUALITY_IDX][:match.start()]]

        # skip reads that are less than 36 bases after trimming
        if len(trimmed_read[SEQUENCE_IDX]) < 36:
            continue

        # Otherwise add trimmed read to output
        reads_without_feature.write("\n".join(trimmed_read) + "\n")

        # Also save complete (untrimmed) reads containing the matched sequence.
        # By mapping these reads to the genome we can eliminate false hits; 
        # i.e. reads that contain a portion of the sequence of intereste but
        # are not actual trans-splicing / poly-adenylation reads.
        untrimmed_read = [read[ID_IDX],
                          read[SEQUENCE_IDX],
                          "+",
                          read[QUALITY_IDX]]
        reads_with_feature.write("\n".join(untrimmed_read) + "\n")

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
    gzip_str(output_without_feature, reads_without_feature)
    gzip_str(output_with_feature, reads_with_feature)
    gzip_str(output_mated_reads, mated_reads_buffer)

    # clean up
    fastq.close()
    fastq_mated.close()
    reads_without_feature.close()
    reads_with_feature.close()
    mated_reads_buffer.close()

    log.info("# Finished processing %s" % os.path.basename(input_file))

def remove_false_hits(feature_name, build_dir, sample_id, read_num):
    """Remove reads that map to genome before trimming"""
    output_dir = ('%s/%s/tophat/%s_unfiltered_untrimmed' % (
        build_dir, sample_id, read_num))
    genome = os.path.splitext(args.genome)[0]

    # input read base directory
    basedir = '%s/%s/fastq' % (build_dir, sample_id)

    # R1 filepath (including matched sequence)
    if (read_num == 'R1'):
        r1_filepath = ('%s/unfiltered/%s_R1_1_with_%s.fastq.gz' %
                       (basedir, sample_id, feature_name))
        r2_filepath = r1_filepath.replace('1_with_%s' % feature_name, '2')
    # R2
    else:
        r2_filepath = ('%s/unfiltered/%s_R2_2_with_%s.fastq.gz' %
                       (basedir, sample_id, feature_name))
        r1_filepath = r2_filepath.replace('2_with_%s' % feature_name, '1')

    # Map reads using Tophat
    # @TODO parameterize extra_args (except for --no-mixed in this case) to
    # allow for easier customization
    loggers[sample_id][feature_name][read_num].info(
        "# Mapping full reads containing sequence match to find false\n"
        "# hits (reads that correspond to actual features in the genome)"
    )
    ret = run_tophat(output_dir, genome, 
                     loggers[sample_id][feature_name][read_num],
                     r1_filepath, r2_filepath,
                     extra_args='--mate-inner-dist 170 --no-mixed')

    # Make sure tophat succeeded
    if ret != 0:
        logging.error(
            "# Error running tophat 1/2! %s (%s)" % (sample_id, read_num)
        )
        sys.exit()

    # Get ids of actual feature-containing reads (those that failed to map when
    # the putative SL/Poly(A) sequence was included).
    sam = pysam.Samfile(os.path.join(output_dir, 'unmapped.bam'), 'rb')
    good_ids = [x.qname for x in sam]

    # number of reads before filtering
    num_reads_before = num_lines(r1_filepath) / 4
    loggers[sample_id][feature_name][read_num].info(
        "# Removing %d false hits (%d total)" %
        (num_reads_before - len(good_ids), num_reads_before)
    )

    # Create true hits directory
    hits_dir = os.path.join(basedir, 'filtered') 

    # Create filtered versions of R1 (and R2) fastq files with only the un-
    # mapped reads
    loggers[sample_id][feature_name][read_num].info(
        "# Filtering matched reads to remove false hits"
    )

    # Filepaths
    r1_infile = r1_filepath.replace('with', 'without')
    r2_infile = r2_filepath.replace('with', 'without')
    r1_outfile = r1_infile.replace('unfiltered', 'filtered')
    r2_outfile = r2_infile.replace('unfiltered', 'filtered')

    filter_fastq(r1_infile, r2_infile, r1_outfile, r2_outfile,
                 good_ids, loggers[sample_id][feature_name][read_num])

    loggers[sample_id][feature_name][read_num].info("# Finished removing false hits.")

def map_reads(feature_name, build_dir, sample_id, read_num):
    """Maps the filtered reads back to the genome"""
    output_dir = '%s/%s/tophat/%s_filtered_trimmed' % (
        build_dir, sample_id, read_num
    )
    genome = os.path.splitext(args.genome)[0]

    loggers[sample_id][feature_name][read_num].info(
        "# Mapping filtered reads back to genome"
    )

    # input read base directory
    basedir = '%s/%s/fastq' % (build_dir, sample_id)

    # R1 filepath (including matched sequence)
    if (read_num == 'R1'):
        r1_filepath = (
            '%s/filtered/%s_R1_1_without_%s.fastq.gz' %
            (basedir, sample_id, feature_name)
        )

        # R2 filepath (for PE reads)
        r2_filepath = r1_filepath.replace('1_with_%s' % feature_name, '2')

        # If SE, set filepath to empty string
        if not os.path.exists(r2_filepath):
            r2_filepath = ""
    # R2
    else:
        r2_filepath = (
            '%s/filtered/%s_R2_2_without_%s.fastq.gz' %
            (basedir, sample_id, feature_name)
        )
        r1_filepath = r2_filepath.replace('2_with_%s' % feature_name, '1')

    # Map reads using Tophat
    #  --no-mixed ?
    ret = run_tophat(output_dir, genome,
                 loggers[sample_id][feature_name][read_num],
                 r1_filepath, r2_filepath,
                 extra_args='--mate-inner-dist 170 --transcriptome-max-hits 1')

    # Make sure tophat succeeded
    if ret != 0:
        logging.error(
            "# Error running tophat 2/2! %s (%s)" % (sample_id, read_num)
        )
        sys.exit()

    loggers[sample_id][feature_name][read_num].info(
        "# Finished mapping hits to genome"
    )


def compute_coordinates(feature_name, build_dir, sample_id, read_num):
    """
    Outputs coordinates and frequencies of putative UTR features to a GFF file.

    References
    ----------
    * http://www.cgat.org/~andreas/documentation/pysam/api.html#pysam.Samfile
    * http://biopython.org/wiki/GFF_Parsing
    * http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html
    """
    # Load existing gene annotations
    annotations_fp = open(args.gff)

    # Get chromosomes from GFF file
    chromosomes = {}

    # output directory
    output_dir = '%s/%s/results' % (build_dir, sample_id)

    for entry in GFF.parse(annotations_fp):
        if len(entry.features) > 0 and entry.features[0].type == 'chromosome':
            chromosomes[entry.id] = entry

    # Create a dictionary to keep track of the coordinates.
    results = {}

    # Create output CSV writer for hits that are not near any genes
    no_nearby_genes = csv.writer(
        open('%s/no_nearby_genes_%s.csv' % (output_dir, read_num), 'w')
    )
    no_nearby_genes.writerow(['read_id', 'chromosome', 'strand', 'position'])

    # Create output CSV writer for hits that are found inside known CDS's
    inside_cds = csv.writer(
        open('%s/inside_cds_%s.csv' % (output_dir, read_num), 'w')
    )
    inside_cds.writerow(['read_id', 'chromosome', 'strand', 'position'])

    # Bam input
    input_bam = '%s/%s/tophat/%s_filtered_trimmed/accepted_hits_sorted.bam' % (
        build_dir, sample_id, read_num
    )

    # file to save results for individual sample
    sample_csv_writer = csv.writer(
        open('%s/matched_reads_%s.csv' % (output_dir, read_num), 'w')
    )
    sample_csv_writer.writerow(['read_id', 'gene_id', 'chromosome',
                                'strand', 'position', 'distance'])

    # keep track of how many reads were found in the expected location
    num_good = 0
    num_no_nearby_genes = 0
    num_inside_cds = 0

    # open sam file
    sam = pysam.Samfile(input_bam, 'rb')

    # Keep track of read id so we only count each one once
    read_ids = []

    # Get coordinate and strand for each read in bam file
    for read in sam:
        # if we have already counted the read, stop here
        if read.qname in read_ids:
            continue
        read_ids.append(read.qname)

        # Chromosome and strand where the read was mapped
        chromosome = sam.getrname(read.tid)
        strand = "-" if read.is_reverse else "+"

        # First, check to make sure the acceptor site does not fall within
        # a known CDS: if it does, save to a separate file to look at later
        if is_inside_cds(chromosomes, read.pos):
            inside_cds.writerow([read.qname, chromosome, strand, read.pos])
            num_inside_cds = num_inside_cds + 1
            continue

        # Find nearest gene
        gene = find_closest_gene(chromosomes, strand, read.pos)

        # If no nearby genes were found, stop here
        if gene is None:
            no_nearby_genes.writerow(
                [read.qname, chromosome, strand, read.pos])
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
            read.qname, gene['id'], chromosome, strand, read.pos,
            gene['distance']
        ])

        # Increment site count and save distance from gene
        if not read.pos in results[chromosome][gene['id']]:
            results[chromosome][gene['id']][read.pos] = {
                "count": 1,
                "distance": gene['distance'],
                "strand": strand,
                "description": gene_description
            }
        else:
            results[chromosome][gene['id']][read.pos]['count'] += 1

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
        output_dir, feature_nam, read_num
    )

    feature_type = 'trans_splice_site' if feature_name is 'sl' else 'polyA_site'
    output_coordinates(results, feature_type, output_filepath)

    logging.info("# Finished!")

    # clean up
    annotations_fp.close()

def is_inside_cds(chromosomes, location):
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
    """
    # Number of bases before or after feature
    half_window = args.window_size / 2

    # Extract region of sequence surrounding the putative feature of
    # interest
    nearby = chromosomes[chromosome][location - half_window:location + half_window]

    # scan all genes near putative feature
    for feature in nearby.features:
        # (half_window is the center point in new sub-region)
        if ((feature.location.start <= half_window) and 
            (feature.location.end >= half_window)):
            return True

    return False

def find_closest_gene(chromosomes, strand, location):
    """
    Finds the closest gene to a specified location that is in the expected
    orientation.
    """
    # 1. Get genes within +/- N bases of location (if any)
    # 2. Find closest match
    if strand == "+":
        # For positive-strand sites, search region just downstream 
        # of splice-site for genes
        subseq = chromosomes[chromosome][location:location + args.window_size]
        gene_start = 'start'
        offset = 0
    else:
        # For negative-strand sites, search region just upstream
        # of splice-site for genes
        subseq = chromosomes[chromosome][location - args.window_size:location]
        gene_start = 'end'
        offset = args.window_size

    # If there are no nearby genes, stop here
    if len(subseq.features) == 0:
        return None

    # Otherwise find closest gene to the acceptor site
    closest_gene = None
    closest_index = None
    closest_dist = float('inf')

    for i, gene in enumerate(subseq.features):
        dist = abs(offset - int(getattr(gene.location, gene_start)))
        if dist < closest_dist:
            closest_index = i
            closest_gene = gene.id
            closest_dist = dist

    # Get description of closest matching gene
    gene_description = (subseq.features[closest_index]
                                .qualifiers['description'].pop())

    # Return details for matching gene
    return {
        'id': closest_gene,
        'description': gene_description,
        'distance': closest_dist
    }

def output_coordinates(results, feature_type, filepath):
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
    fp.write("##gff-version\t3")
    fp.write("##feature-ontology\tsofa.obo")
    fp.write("##attribute-ontology\tgff3_attributes.obo")

    # Write header to output
    writer = csv.writer(fp, delimiter='\t')

    # write output to csv
    for chrnum in results:
        for gene_id in results[chrnum]:
            for acceptor_site in results[chrnum][gene_id]:
                # gff3 attributes
                attributes = "ID=%s;Name=%s;description=%s" % (
                    gene_id, gene_id,
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

# Poly(A) tail sub-directory
polya_build_dir = os.path.join(
    args.build_directory,
    'poly-a',
    'minlength-%d' % args.min_polya_length,
    'anchored' if args.exclude_internal_matches else 'unanchored'
)

# Poly(A) tail reverse complement sub-directory
polyt_build_dir = os.path.join(
    args.build_directory,
    'poly-t',
    'minlength-%d' % args.min_polya_length,
    'anchored' if args.exclude_internal_matches else 'unanchored'
)

# Get a list of sample ids ids
# @TODO: Generalize handling of sample ids
input_regex = re.compile(r'.*(HPGL[0-9]+).*')
sample_ids = []

# currently processing
for filename in glob.glob(args.input_reads):
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
for sample_id in sample_ids:
    # create build directories
    for base_dir in [sl_build_dir, polya_build_dir, polyt_build_dir]:
        for sub_dir in ['fastq/filtered', 'fastq/unfiltered', 'ruffus', 
                        'results', 'log', 'tophat']:
            outdir = os.path.join(base_dir, sample_id, sub_dir)
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
for sample_id in sample_ids_all:
    loggers[sample_id] = {}

    for analysis in ['sl', 'polya', 'polyt']:
        loggers[sample_id][analysis] = {}

        for read_num in ['R1', 'R2']:
            if analysis == 'sl':
                bdir = sl_build_dir
            elif analysis == 'polya':
                bdir = polya_build_dir
            else:
                bdir = polyt_build_dir

            sample_log_name = get_next_log_name(
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
# Ruffus tasks
#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
# Step 1: Find reads with sequence of interest
#
# Finds reads containing a minimum number of bases of the feature of interest;
# either a portion of the spliced leader (SL) sequence, or a Poly(A) or Poly(T)
# tract in the expected location for a polyadenylation event.
#-----------------------------------------------------------------------------

# Input regular expression
#
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
# \1 - directory
# \2 - HPGLxxxx
# \3 - _anything_between_id_and_read_num_
# \4 - R1/R2
# \5 - _anything_after_read_num_
@follows(check_for_bowtie_index)
@follows(check_for_genome_fasta)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq'),
           r'%s/\2/ruffus/\2_\4.find_sl_reads' % sl_build_dir,
           r'\2', r'\4')
def find_sl_reads(input_file, output_file, sample_id, read_num):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing a subsequence of the spliced leader.
    """
    # First-stage filter string
    sl_filter = args.spliced_leader[-args.min_sl_length:]

    # Determine strings to match in reads
    if args.exclude_internal_matches:
        sl_regex = '|'.join(
            ["^" + args.spliced_leader[-x:] for x in 
             range(args.min_sl_length, len(args.spliced_leader) + 1)]
        )
    else:
        sl_regex = sl_filter

    find_sequence(input_file, 'sl', sl_filter, sl_regex, sl_build_dir, 
                  sample_id, read_num)

    # Let Ruffus know we are done
    open(output_file, 'w').close()

@follows(find_sl_reads)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq'),
           r'%s/\2/ruffus/\2_\4.find_polya_reads' % polya_build_dir,
           r'\2', r'\4')
def find_polya_reads(input_file, output_file, sample_id, read_num):
    """Matches reads with possible Poly(A) tail fragment"""
    # Match reads with at least n A's at the end of the read; For now we will 
    # always require matches to be at the end of the read.
    polya_filter = 'A' * args.min_polya_length
    polya_regex = 'A{%d,}$' % (args.min_polya_length)

    find_sequence(input_file, 'polya', polya_filter, polya_regex, 
                  polya_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

@follows(find_polya_reads)
@transform(args.input_reads,
           regex(r'^(.*/)?(HPGL[0-9]+)_(.*)(R[1-2])_(.+)\.fastq'),
           r'%s/\2/ruffus/\2_\4.find_polyt_reads' % polyt_build_dir,
           r'\2', r'\4')
def find_polyt_reads(input_file, output_file, sample_id, read_num):
    """Matches reads with possible Poly(T) tail fragment"""
    # Match reads with at least n T's at the beginning of the read; For now 
    # we will always require matches to be at the beginning of the read.
    polyt_filter = 'T' * args.min_polya_length
    polyt_regex = 'T{%d,}$' % (args.min_polya_length)
    find_sequence(input_file, 'polyt', polyt_filter, polyt_regex, 
                  polyt_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Step 2: Remove false hits
#
# In this step, the matching reads from step 1 are mapped as-is: any reads
# which successfully map in this way can then be attributed to sequences found
# in the actual genome, and not a trans-splicing or polyadenylation event, and
# are filtered out.
#-----------------------------------------------------------------------------
@transform(find_polyt_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).find_polyt_reads'),
           r'\1/\2_\3.remove_sl_false_hits',
           r'\2', r'\3')
def remove_sl_false_hits(input_file, output_file, sample_id, read_num):
    remove_false_hits('sl', sl_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

@transform(remove_sl_false_hits,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).remove_sl_false_hits'),
           r'\1/\2_\3.remove_polya_false_hits',
           r'\2', r'\3')
def remove_polya_false_hits(input_file, output_file, sample_id, read_num):
    """Remove Poly(A) reads that map to genome before trimming"""
    remove_false_hits('polya', polya_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

@transform(remove_polya_false_hits,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).remove_polya_false_hits'),
           r'\1/\2_\3.remove_polyt_false_hits',
           r'\2', r'\3')
def remove_polyt_false_hits(input_file, output_file, sample_id, read_num):
    """Remove Poly(T) reads that map to genome before trimming"""
    remove_false_hits('polyt', polyt_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Step 3: Map trimmed reads
#
# Trims the matched sequences from reads and map to genome. For reads where
# the matched sequence comes from a trans-splicing or polyadenylation event,
# the location of the mapped trimmed read is where the addition took place.
#-----------------------------------------------------------------------------
@transform(remove_polyt_false_hits,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).remove_polyt_false_hits'),
           r'\1/\2_\3.map_sl_reads',
           r'\2', r'\3')
def map_sl_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered spliced-leader containing reads back to the genome"""
    map_reads('sl', sl_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

@transform(map_sl_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).map_sl_reads'),
           r'\1/\2_\3.map_polya_reads',
           r'\2', r'\3')
def map_polya_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered poly-adenylated reads back to the genome"""
    map_reads('polya', polya_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

@transform(map_polya_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).map_polya_reads'),
           r'\1/\2_\3.map_polyt_reads',
           r'\2', r'\3')
def map_polyt_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered Poly(T) reads back to the genome"""
    map_reads('polyt', polyt_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Step 4: Compute UTR coordinates
#
# In this step, the locations of mapped trimmed reads above are used to
# determine likely trans-splicing and polyadenylation acceptor sites in the
# genome for the processed input samples. A couple of additional filtering
# steps are based to ensure that the putative acceptor sites lie in acceptable
# locations (e.g. not inside a CDS.)
#-----------------------------------------------------------------------------

#@collate(map_polyt_reads,
#       regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).map_polyt_reads'),
#       r'\1/\2_.compute_sl_coordinates')
@transform(map_polyt_reads,
           regex(r'^(.*)/(HPGL[0-9]+)_(R[12]).map_polyt_reads'),
           r'\1/\2_\3.compute_sl_coordinates',
           r'\2', r'\3')
def compute_sl_coordinates(input_file, output_file, sample_id, read_num):
    logging.info("# Computing coordinates for mapped trans-splicing events")
    compute_coordinates('sl', sl_build_dir, sample_id, read_num)
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Run pipeline
#-----------------------------------------------------------------------------
if __name__ == "__main__":
    pipeline_run([compute_sl_coordinates], logger=logging.getLogger(''),
                 multiprocess=args.num_threads)
    pipeline_printout_graph("utr_analysis_flowchart.png", "png",
                            [compute_sl_coordinates])

