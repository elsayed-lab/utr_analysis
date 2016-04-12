#!/bin/env python2
# -*- coding: utf-8 -*-
"""
Trypanosomatid Untranslated Region (UTR) analysis Pipeline
Keith Hughitt (khughitt@umd.edu)

Overview
--------
The purpose of this script is the scan a collection of RNA-Seq reads and
determine the location of the 5'UTR splice acceptor sites or 3'UTR
polyadenylation sites.

A reference genome and CDS coordinates are required as input and will be used
to map reads back after removing the spliced leader in order to determine the
SL acceptor site.
"""
import os
import re
import glob
import time
import logging
import warnings
import subprocess
from ruffus import *
from Bio import Seq,BiopythonDeprecationWarning
from utr.io import parse_input, create_build_dirs
from utr.mapping import map_reads, filter_mapped_reads
from utr.util import get_next_file, setup_loggers
from utr.sequence import find_sequence, compute_coordinates

# Hide Biopython deprecation warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

#--------------------------------------
# Global variables
#--------------------------------------

# parse input
args = parse_input()
build_dirs = create_build_dirs(args, sample_ids)

# Get a list of sample ids, e.g. /path/to/input/samples/sample01/...

# For simplicity, the sample ID for a given file is assumed to be everything
# in the filename up to the mate pair specification.
input_regex = re.compile(r'^(.*/)?(.*)_(R?[1-2])')

SAMPLE_REGEX_DIR_IDX = 0
SAMPLE_REGEX_ID_IDX = 1
SAMPLE_REGEX_READ_NUM_IDX = 2

sample_ids = []

for filename in glob.glob(args.input_reads):
    sample_id = re.match(input_regex, filename).groups()[SAMPLE_REGEX_ID_IDX]
    if sample_id not in sample_ids:
        sample_ids.append(sample_id)

# setup global and task-specific loggers
loggers = setup_loggers(args.build_directory, sample_ids)

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
# Ex. "$RAW/tcruzir21/Sample0121/processed/Sample0121_R1_filtered.fastq"
#
#    \1 - directory                   ($RAW/tcruzir21/Sample0121/processed/)
#    \2 - sample id / filename prefix (Sample0121)
#    \3 - read number                 (R1)
#    \4 - _anything_after_read_num_   (filtered)
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
           regex(r'^(.*/)?([^_]*)_(R?[1-2])_?(.*)?\.fastq(\.gz)?'),
           r'%s/\2/ruffus/\2_\3.filter_nontarget_reads' % build_dirs['shared'],
           r'\2', r'\3')
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
    tophat_dir = os.path.join(build_dirs['shared'], sample_id,
                          'tophat', 'mapped_to_nontarget')
    fastq_dir = os.path.join(build_dirs['shared'], sample_id,
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
# The next step is to remove reads which map to the target genome before any
# trimming. These reads likely represent actual features in the genome and not
# the capped or polyadenylated reads we are interested in.
#-----------------------------------------------------------------------------
#@transform(args.input_reads,
#           regex(r'^(.*/)?([^_]*)_(R?[1-2])_?(.*)?\.fastq(\.gz)?'),
#           r'%s/\2/ruffus/\2_\3.filter_genomic_reads' % build_dirs['shared'],
#           r'\2', r'\3')
@follows(filter_nontarget_reads)
@subdivide(args.input_reads,
           regex(r'^(.*/)?([^_]*)_(R?[1-2])_?(.*)?\.fastq(\.gz)?'),
           [r'%s/\2/ruffus/\2_\3.filter_genomic_reads.sl' %
               build_dirs['shared'],
            r'%s/\2/ruffus/\2_\3.filter_genomic_reads.rsl' %
            build_dirs['shared'],
            r'%s/\2/ruffus/\2_\3.filter_genomic_reads.polya' %
            build_dirs['shared'],
            r'%s/\2/ruffus/\2_\3.filter_genomic_reads.polyt' %
            build_dirs['shared']],
            r'%s/\2/ruffus/\2_\3.filter_genomic_reads' % build_dirs['shared'],
            r'\2', r'\3')
def filter_genomic_reads(input_file, output_files, output_base, sample_id, read_num):
    # Clean-up any files from previous runs
    for x in output_files:
        if os.path.exists(x):
            os.unlink(x)

    # We only need to map once for each mated pair
    if read_num == "R2":
        # Wait for R1 task to finish processing and then mark as finished
        while not os.path.exists(output_files[0].replace("R2", "R1")):
            time.sleep(120)

        # Mark as finished and exit
        for output_file in output_files:
            open(output_file, 'w').close()
        return

    logging.info("# Removing genomic reads.")

    # determine source of input reads to use
    if args.nontarget_genome:
        # if nontarget genome was specified, use output from that step
        input_fastq_dir = os.path.join(build_dirs['shared'], sample_id,
                                    'fastq', 'nontarget_reads_removed')
        r1 = os.path.join(input_fastq_dir,
                        "%s_nontarget_reads_removed.1.fastq.gz" % (sample_id))
        r2 = r1.replace(".1", ".2")
    else:
        # otherwise use inputs specified at run-time
        r1 = input_file
        r2 = r1.replace("R1", "R2")

    # output tophat and fastq directories
    tophat_dir = os.path.join(build_dirs['shared'], sample_id,
                              'tophat', 'mapped_to_target_untrimmed')
    output_fastq_dir = os.path.join(build_dirs['shared'], sample_id,
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
    for output_file in output_files:
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
           regex(r'^(.*)/(.*)_(R?[12]).filter_genomic_reads.sl'),
           r'%s/\2/ruffus/\2_\3.find_sl_reads' % build_dirs['sl'],
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
        build_dirs['shared'], sample_id, 'fastq', 'genomic_reads_removed',
        "%s_genomic_reads_removed.%s.fastq.gz" % (sample_id, read_num[-1]))

    log_handle = loggers[sample_id]['sl'][read_num]

    find_sequence(input_reads, 'sl', sl_filter, sl_regex, build_dirs['sl'],
                  sample_id, read_num, log_handle)

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
           regex(r'^(.*)/(.*)_(R?[12]).find_sl_reads'),
           r'\1/\2_\3.map_sl_reads',
           r'\2', r'\3')
def map_sl_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered spliced-leader containing reads back to the genome"""
    log_handle = loggers[sample_id]['sl'][read_num]
    map_reads('sl', build_dirs['sl'], sample_id, read_num, log_handle)
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
           regex(r'^(.*)/(.*)_(R?[12]).map_sl_reads'),
           r'\1/\2_\3.compute_sl_coordinates',
           r'\2', r'\3')
def compute_sl_coordinates(input_file, output_file, sample_id, read_num):
    logging.info("# Computing coordinates for mapped trans-splicing events "
                 "(%s)" % sample_id)
    compute_coordinates('sl', build_dirs['sl'], sample_id, read_num)
    logging.info("# Finished!")    
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Reverse SL analysis
#-----------------------------------------------------------------------------

#
# Reverse SL Step 1
#
@transform(filter_genomic_reads,
           regex(r'^(.*)/(.*)_(R?[12]).filter_genomic_reads.rsl'),
           r'%s/\2/ruffus/\2_\3.find_rsl_reads' % build_dirs['rsl'],
           r'\2', r'\3')
def find_rsl_reads(input_file, output_file, sample_id, read_num):
    """Matches reads with possible Poly(A) tail fragment"""
    # Create a reversed version of the SL sequence
    reverse_sl = str(Seq.Seq(args.spliced_leader).reverse_complement())

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
        build_dirs['shared'], sample_id, 'fastq', 'genomic_reads_removed',
        "%s_genomic_reads_removed.%s.fastq.gz" % (sample_id, read_num[-1]))

    logging.info("# find_rsl_reads (%s %s)" % (sample_id, read_num))

    find_sequence(input_reads, 'rsl', rsl_filter, rsl_regex,
                  build_dirs['rsl'], sample_id, read_num, trimmed_side='right',
                  reverse=True)
    open(output_file, 'w').close()

#
# RSL Step 2
#
@transform(find_rsl_reads,
           regex(r'^(.*)/(.*)_(R?[12]).find_rsl_reads'),
           r'\1/\2_\3.map_rsl_reads',
           r'\2', r'\3')
def map_rsl_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered poly-adenylated reads back to the genome"""
    log_handle = loggers[sample_id]['rsl'][read_num]
    map_reads('rsl', build_dirs['rsl'], sample_id, read_num, log_handle)
    open(output_file, 'w').close()

#
# RSL Step 3
#
@transform(map_rsl_reads,
           regex(r'^(.*)/(.*)_(R?[12]).map_rsl_reads'),
           r'\1/\2_\3.compute_rsl_coordinates',
           r'\2', r'\3')
def compute_rsl_coordinates(input_file, output_file, sample_id, read_num):
    logging.info(
        "# Computing coordinates for mapped reverse sl events [A]")
    compute_coordinates('rsl', build_dirs['rsl'], sample_id, read_num)
    logging.info("# Finished!")    
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Poly(A) Analysis
#-----------------------------------------------------------------------------

#
# Poly(A) Step 1
#
@transform(filter_genomic_reads,
           regex(r'^(.*)/(.*)_(R?[12]).filter_genomic_reads.polya'),
           r'%s/\2/ruffus/\2_\3.find_polya_reads' % build_dirs['polya'],
           r'\2', r'\3')
def find_polya_reads(input_file, output_file, sample_id, read_num):
    """Matches reads with possible Poly(A) tail fragment"""
    polya_filter = 'A' * args.min_polya_length

    if args.exclude_internal_polya_matches:
        polya_regex = 'A{%d,}$' % (args.min_polya_length)
    else:
        polya_regex = 'A{%d,}' % (args.min_polya_length)

    # input reads
    input_reads = os.path.join(
        build_dirs['shared'], sample_id, 'fastq', 'genomic_reads_removed',
        "%s_genomic_reads_removed.%s.fastq.gz" % (sample_id, read_num[-1]))

    logging.info("# find_polya_reads (%s %s)" % (sample_id, read_num))

    find_sequence(input_reads, 'polya', polya_filter, polya_regex,
                  build_dirs['polya'], sample_id, read_num, trimmed_side='right')
    open(output_file, 'w').close()

#
# Poly(A) Step 2
#
@transform(find_polya_reads,
           regex(r'^(.*)/(.*)_(R?[12]).find_polya_reads'),
           r'\1/\2_\3.map_polya_reads',
           r'\2', r'\3')
def map_polya_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered poly-adenylated reads back to the genome"""
    log_handle = loggers[sample_id]['polya'][read_num]
    map_reads('polya', build_dirs['polya'], sample_id, read_num, log_handle)
    open(output_file, 'w').close()

#
# Poly(A) Step 3
#
@transform(map_polya_reads,
           regex(r'^(.*)/(.*)_(R?[12]).map_polya_reads'),
           r'\1/\2_\3.compute_polya_coordinates',
           r'\2', r'\3')
def compute_polya_coordinates(input_file, output_file, sample_id, read_num):
    logging.info(
        "# Computing coordinates for mapped polyadenylation events [A]")
    compute_coordinates('polya', build_dirs['polya'], sample_id, read_num)
    logging.info("# Finished!")    
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Poly(T) Analysis
#-----------------------------------------------------------------------------

#
# Poly(T) Step 1
#
@transform(filter_genomic_reads,
           regex(r'^(.*)/(.*)_(R?[12]).filter_genomic_reads.polyt'),
           r'%s/\2/ruffus/\2_\3.find_polyt_reads' % build_dirs['polyt'],
           r'\2', r'\3')
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
        build_dirs['shared'], sample_id, 'fastq', 'genomic_reads_removed',
        "%s_genomic_reads_removed.%s.fastq.gz" % (sample_id, read_num[-1]))

    logging.info("# find_polyt_reads (%s %s)" % (sample_id, read_num))

    find_sequence(input_reads, 'polyt', polyt_filter, polyt_regex,
                  build_dirs['polyt'], sample_id, read_num, reverse=True)
    open(output_file, 'w').close()

#
# Poly(T) Step 2
#
@transform(find_polyt_reads,
           regex(r'^(.*)/(.*)_(R?[12]).find_polyt_reads'),
           r'\1/\2_\3.map_polyt_reads',
           r'\2', r'\3')
def map_polyt_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered Poly(T) reads back to the genome"""
    log_handle = loggers[sample_id]['polyt'][read_num]
    map_reads('polyt', build_dirs['polyt'], sample_id, read_num, log_handle)
    open(output_file, 'w').close()

#
# Poly(T) Step 3
#
@transform(map_polyt_reads,
           regex(r'^(.*)/(.*)_(R?[12]).map_polyt_reads'),
           r'\1/\2_\3.compute_polyt_coordinates',
           r'\2', r'\3')
def compute_polyt_coordinates(input_file, output_file, sample_id, read_num):
    logging.info("# Computing coordinates for mapped polyadenylation events [T]")
    log_handle = loggers[sample_id]['polyt'][read_num]
    compute_coordinates('polyt', build_dirs['polyt'], sample_id, read_num,
                        log_handle)
    logging.info("# Finished!")    
    open(output_file, 'w').close()

#-----------------------------------------------------------------------------
# Step 6: Combine coordinate output for multiple samples
#
# Next, we create a single GFF file using the information contained in the GFF
# files that were constructed for each sample.
#
# Coordinates with low coverage may also be filtered out at this step.
#-----------------------------------------------------------------------------
@merge(compute_polyt_coordinates, '%s/finished' % build_dirs['combined'])
def combine_results(input_files, output_file):
    # Convert input ruffus tasks to corresponding GFF filepaths
    regex = '.*/(.*)_(R?[1-2]).*'

    # Combine spliced leader output
    logging.info("# Combining spliced leader coordinates output")
    sl_gffs = []
    for infile in input_files:
        (sample_id, read_num) = re.match(regex, infile).groups()
        # SL
        gff1 = "%s/%s/results/sl_coordinates_%s.gff" % (
                 build_dirs['sl'], sample_id, read_num)
        sl_gffs.append(gff1)

        # Reverse SL
        gff2 = "%s/%s/results/rsl_coordinates_%s.gff" % (
                 build_dirs['rsl'], sample_id, read_num)
        sl_gffs.append(gff2)

    sl_outfile = os.path.join(build_dirs['combined'], 'spliced_leader.gff')
    sl_combined = combine_gff_results(sl_gffs)

    # Combine Poly(A) output
    logging.info("# Combining Poly(A) coordinates output")

    polya_gffs = []
    for infile in input_files:
        (sample_id, read_num) = re.match(regex, infile).groups()
        # Poly(A)
        gff1 = "%s/%s/results/polya_coordinates_%s.gff" % (
                 build_dirs['polya'], sample_id, read_num)
        polya_gffs.append(gff1)

        # Poly(T)
        gff2 = "%s/%s/results/polyt_coordinates_%s.gff" % (
                 build_dirs['polyt'], sample_id, read_num)
        polya_gffs.append(gff2)

    polya_outfile = os.path.join(build_dirs['combined'], 'polya.gff')
    polya_combined = combine_gff_results(polya_gffs)

    # Save summary GFFs
    output_coordinates(sl_combined, 'sl', sl_outfile, track_color='83,166,156')
    output_coordinates(polya_combined, 'polya', polya_outfile,
                       track_color='166,83,93')

#-----------------------------------------------------------------------------
# Run pipeline
#-----------------------------------------------------------------------------
if __name__ == "__main__":
    #pipeline_run([compute_rsl_coordinates], forcedtorun_tasks=[combine_results],
    #             touch_files_only=True, logger=logging.getLogger(''))
    pipeline_run([combine_results], forcedtorun_tasks=[combine_results],
                 logger=logging.getLogger(''), multiprocess=args.num_threads)
    pipeline_printout_graph("utr_analysis_flowchart.png", "png",
                            [combine_results])

