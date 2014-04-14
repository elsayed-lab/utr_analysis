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
import StringIO
import textwrap
import subprocess
from ruffus import *
from Bio import Seq,SeqIO
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
                        help='RNA-Seq FASTQ or gzipped FASTQ glob string')
    parser.add_argument('-d', '--build-directory', required=True,
                        help='Directory to save output to')
    parser.add_argument('-f1', '--target-genome', dest='target_genome',
                        required=True, help=('Genome sequence FASTA filepath '
                        'for target species'))
    parser.add_argument('-f2', '--nontarget-genome', dest='nontarget_genome', 
                        help=('Genome sequence FASTA filepath for species to '
                              'be filtered out prior to mapping. (optional)'))
    parser.add_argument('-g', '--gff-annotation', dest='gff', required=True,
                        help='Genome annotation GFF')
    parser.add_argument('-s', '--sl-sequence', dest='spliced_leader',
                        required=True, help='Spliced leader DNA sequence', 
                        default=None)
    parser.add_argument('-es', '--exclude-internal-sl-matches',
                        help=('Only allow matches with the SL at the upstream'
                              'end of a read.'), action='store_true')
    parser.add_argument('-ep', '--exclude-internal-polya-matches',
                        help=('Only allow matches with the Poly(A) tail at the'
                              'downstream end of a read.'), action='store_true')
    parser.add_argument('-m', '--min-sl-length', default=10, type=int,
                        help='Minimum length of SL match (default=10)')
    parser.add_argument('-p', '--min-polya-length', default=10, type=int,
                        help='Minimum length of Poly-A match (default=10)')
    parser.add_argument('-w', '--window-size', default=10000, type=int,
                        help=('Number of bases up or downstream of feature to'
                              'scan for related genes (default=10000)'))
    parser.add_argument('--num-threads', default=4, type=int,
                        help='Number of threads to use.')

    # Parse arguments
    args = parser.parse_args()

    # Replace any environmental variables and return args
    args.target_genome = os.path.expandvars(args.target_genome)
    args.input_reads = os.path.expandvars(args.input_reads)
    args.gff = os.path.expandvars(args.gff)

    if args.nontarget_genome:
        args.nontarget_genome = os.path.expandvars(args.nontarget_genome)

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

    # FASTA filename without extension
    genome_basename = os.path.splitext(genome)[0]

    # @NOTE 2014/02/06 -- Tophat 2.0.10 fails for some reads when attempting
    # to use more than one thread. For now, just perform the mapping
    # using a single-thread to be safe...

    # build command
    cmd = "tophat --num-threads %d --max-multihits %d %s -o %s %s %s %s" % (
           num_threads, max_multihits, extra_args, 
           output_dir, genome_basename, r1, r2)

    # run tophat
    ret = run_command(cmd, log_handle)

    # check to see if tophat succeeded and stop execution otherwise
    if ret != 0:
        log_handle.error("# Error running tophat (%s)!" % genome)
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

def filter_mapped_reads(r1, r2, genome, tophat_dir, output_fastq, log_handle):
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
    """
    # bam / fastq filepaths
    bam_input = os.path.join(tophat_dir, 'unmapped_sorted.bam')

    # check to see if tophat output already exists
    if not os.path.exists(bam_input):
        # map reads to genome using tophat
        log_handle.info("# Mapping against %s" % os.path.basename(genome))
        ret = run_tophat(tophat_dir, genome, log_handle, r1, r2,
                        extra_args='--no-mixed')

        # number of reads before filtering
        num_reads_total = num_lines(r1) / 4
        num_unmapped = num_lines(bam_input) / 2

        log_handle.info(
            "# Ignoring %d reads which mapped to specified genome (%d total)" %
            (num_reads_total - num_unmapped, num_unmapped)
        )
    else:
        log_handle.info("# Skipping %s: output already exists" % tophat_dir)

    # keep unampped reads
    log_handle.info("# Converting Tophat bam output to FASTQ")
    run_bam2fastx(bam_input, output_fastq, log_handle)

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
                  build_dir, sample_id, read_num, trim_direction='left',
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
    trim_direction: str
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

    # list to keep track of potential matches
    matches = []

    # output filepaths
    output_base = '%s/%s/fastq/unfiltered/%s_%s_%s' % (
        build_dir, sample_id, sample_id, read_num, read_num[-1]
    )
    output_untrimmed = "%s_%s_untrimmed.fastq.gz" % (output_base, feature_name)
    output_trimmed = "%s_%s_trimmed.fastq.gz" % (output_base, feature_name)

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
        if sequence_filter not in read[SEQUENCE_IDX]:
            continue

        # check for match
        match = re.search(read_regex, read[SEQUENCE_IDX])

        # move on the next read if no match is found
        if match is None:
            continue

        # otherwise add to output fastq

        # If feature is expected to be found at left of read, trim everything 
        # up to the end of the match
        if trim_direction == 'left':
            trimmed_read = [read[ID_IDX],
                            read[SEQUENCE_IDX][match.end():],
                            "+",
                            read[QUALITY_IDX][match.end():]]
        else:
            # otherwise trim from the start of the match to the end of the read
            trimmed_read = [read[ID_IDX],
                            read[SEQUENCE_IDX][:match.start()],
                            "+",
                            read[QUALITY_IDX][:match.start()]]

        # skip reads that are less than 36 bases after trimming
        if len(trimmed_read[SEQUENCE_IDX]) < 36:
            continue

        # take reverse complement if requested
        if reverse:
            trimmed_read[SEQUENCE_IDX] = str(
                Seq.Seq(trimmed_read[SEQUENCE_IDX]).reverse_complement())
            trimmed_read[QUALITY_IDX] = trimmed_read[QUALITY_IDX][::-1]

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

    log.info("# Finished processing %s" % os.path.basename(input_file))

def map_reads(feature_name, build_dir, sample_id, read_num):
    """Maps the filtered reads back to the genome"""
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
                 r1_filepath, r2_filepath,
                 extra_args='--mate-inner-dist 170 --transcriptome-max-hits 1')

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

    # If processing Poly(A)/Poly(T), load genome sequence as well
    # This will be used to make sure number of A's or T's in read exceeds
    # the number found at the mapped location in the genome.
    if feature_name != 'sl':
        genome = SeqIO.parse(args.target_genome, 'fasta')
        chr_sequences = {x.id:x for x in genome}

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

    # For Poly(A) analysis, we will also load the untrimmed version to
    # determine original read length
    if feature_name in ['polya', 'polyt']:
        input_bam_untrimmed = '%s/%s/tophat/mapped_to_target_untrimmed/unmapped_sorted.bam' % (
            shared_build_dir, sample_id
        )
        sam_untrimmed = pysam.Samfile(input_bam_untrimmed, 'rb')

        read_lengths = {}

        for read in sam_untrimmed:
            rnum = 'R1' if read.is_read1 else 'R2'
            if not read.qname in read_lengths:
                read_lengths[read.qname] = {rnum:read.rlen}
            else:
                read_lengths[read.qname][rnum] = read.rlen

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
        # Get read where feature sequence was found
        #if read.qname in read_ids:
        if ((read.is_read1 and read_num == 'R2') or 
            (read.is_read2 and read_num == 'R1')):
            continue
        read_ids.append(read.qname)

        # Chromosome and strand where the read was mapped
        chromosome = sam.getrname(read.tid)
        strand = "-" if read.is_reverse else "+"

        # SL/Poly(A) acceptor site location
        if ((strand == '+' and feature_name == 'polya') or 
            (strand == '-' and feature_name == 'polyt') or
            (strand == '-' and feature_name == 'sl')):
            # acceptor site at right end of read
            acceptor_site = read.pos + read.rlen
        else:
            # acceptor site at left end of read
            acceptor_site = read.pos

        # First, check to make sure the acceptor site does not fall within
        # a known CDS: if it does, save to a separate file to look at later
        if is_inside_cds(chromosomes[chromosome], acceptor_site):
            inside_cds.writerow([read.qname, chromosome, strand, acceptor_site])
            num_inside_cds = num_inside_cds + 1
            continue

        # For Poly(A)/Poly(T) reads, check to see if reads contain at least 1
        # more A/T at the end of read compared with location mapped in genome
        if feature_name in ['polya', 'polyt']:
            # Number of A's or T's
            feature_length = read_lengths[read.qname][read_num] - read.rlen

            # Count number of A's / T's just downstream of acceptor site
            if ((feature_name == 'polya' and strand == '+') or 
                (feature_name == 'polyt' and strand == '-')):
                # Check for A's at right end of read
                rec = chr_sequences[chromosome][read.pos + read.rlen:read.pos + read.rlen + feature_length]

                # Make sure that read contained at least one more A than is
                # found in the genome at mapped location
                if rec.seq.count('A') >= feature_length:
                    continue
            else:
                # Check for T's at right end of read
                rec = chr_sequences[chromosome][read.pos - feature_length:read.pos]

                # Make sure that read contained at least one more A than is
                # found in the genome at mapped location
                if rec.seq.count('T') >= feature_length:
                    continue

        # Find nearest gene
        gene = find_closest_gene(chromosomes[chromosome], strand, acceptor_site)

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
            read.qname, gene['id'], chromosome, strand, acceptor_site,
            gene['distance']
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

    # clean up
    annotations_fp.close()

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
    """
    # Number of bases before or after feature
    half_window = args.window_size / 2

    # Extract region of sequence surrounding the putative feature of
    # interest
    nearby = chromosome[location - half_window:location + half_window]

    # scan all genes near putative feature
    for feature in nearby.features:
        # (half_window is the center point in new sub-region)
        if ((feature.location.start <= half_window) and 
            (feature.location.end >= half_window)):
            return True

    return False

def find_closest_gene(chromosome, strand, location):
    """
    Finds the closest gene to a specified location that is in the expected
    orientation.
    """
    # 1. Get genes within +/- N bases of location (if any)
    # 2. Find closest match
    if strand == "+":
        # For positive-strand sites, search region just downstream 
        # of splice-site for genes
        subseq = chromosome[location:location + args.window_size]
        gene_start = 'start'
        offset = 0
    else:
        # For negative-strand sites, search region just upstream
        # of splice-site for genes
        subseq = chromosome[location - args.window_size:location]
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
    desc_qualifiers = subseq.features[closest_index].qualifiers['description']
    if len(desc_qualifiers) > 0:
        gene_description = desc_qualifiers[0]
    else:
        # TESTING 2014/03/21 (can remove if no longer occurs)
        print("Missing description for gene: %s" % closest_gene)
        gene_description = ""

    # Return details for matching gene
    return {
        'id': closest_gene,
        'description': gene_description,
        'distance': closest_dist
    }

def output_coordinates(results, feature_name, filepath):
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

    # Write header to output
    writer = csv.writer(fp, delimiter='\t')

    # Determine GFF feature type to use
    feature_type = 'trans_splice_site' if feature_name is 'sl' else 'polyA_site'

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

def combine_gff_results(input_gffs, outfile, feature_name):
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

        # Parse with CSV reader
        reader = csv.DictReader(fp, delimiter='\t', fieldnames=gff_fields)

        for row in reader:
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
            description = parts[-1]

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

    # Save summary GFF
    output_coordinates(results, feature_name, outfile)

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
    'minlength-%d' % args.min_sl_length,
    'anchored' if args.exclude_internal_sl_matches else 'unanchored'
)
if not os.path.exists(combined_output_dir):
    os.makedirs(combined_output_dir, mode=0o755)

# Create a unique build paths for specified parameters

# SL sub-directory
sl_build_dir = os.path.join(
    args.build_directory,
    'spliced_leader',
    'minlength-%d' % args.min_sl_length,
    'anchored' if args.exclude_internal_sl_matches else 'unanchored'
)

# Poly(A) tail sub-directory
polya_build_dir = os.path.join(
    args.build_directory,
    'poly-a',
    'minlength-%d' % args.min_polya_length,
    'anchored' if args.exclude_internal_polya_matches else 'unanchored'
)

# Poly(A) tail reverse complement sub-directory
polyt_build_dir = os.path.join(
    args.build_directory,
    'poly-t',
    'minlength-%d' % args.min_polya_length,
    'anchored' if args.exclude_internal_polya_matches else 'unanchored'
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
    for base_dir in [sl_build_dir, polya_build_dir, polyt_build_dir]:
        for sub_dir in ['fastq/filtered', 'fastq/unfiltered',
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

    # map reads and remove hits
    filter_mapped_reads(r1, r2, args.nontarget_genome, tophat_dir,
                        output_fastq, logging)
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
    filter_mapped_reads(r1, r2, args.target_genome, tophat_dir, output_fastq, 
                        logging)
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

    # Determine strings to match in reads
    if args.exclude_internal_sl_matches:
        sl_regex = '|'.join(
            ["^" + args.spliced_leader[-x:] for x in 
             range(args.min_sl_length, len(args.spliced_leader) + 1)]
        )
    else:
        sl_regex = '|'.join(
            [args.spliced_leader[-x:] for x in 
             range(args.min_sl_length, len(args.spliced_leader) + 1)]
        )

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
# Poly(A) Analysis
#-----------------------------------------------------------------------------

#
# Poly(A) Step 1
#
@follows(compute_sl_coordinates)
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

    find_sequence(input_reads, 'polya', polya_filter, polya_regex, 
                  polya_build_dir, sample_id, read_num, trim_direction='right')
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
        "# Computing coordinates for mapped polyadenylation events [1/2]")
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
        "# Computing coordinates for mapped polyadenylation events [2/2]")
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
        gff = "%s/%s/results/sl_coordinates_%s.gff" % (
                 sl_build_dir, sample_id, read_num)
        sl_gffs.append(gff)

    sl_outfile = os.path.join(combined_output_dir, 'spliced_leader.gff')
    combine_gff_results(sl_gffs, sl_outfile, 'sl')

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
    combine_gff_results(polya_gffs, polya_outfile, 'polya')


#-----------------------------------------------------------------------------
# Run pipeline
#-----------------------------------------------------------------------------
if __name__ == "__main__":
    pipeline_run([combine_results], forcedtorun_tasks=[combine_results],
                 logger=logging.getLogger(''), multiprocess=args.num_threads)
    pipeline_printout_graph("utr_analysis_flowchart.png", "png",
                            [combine_results])

