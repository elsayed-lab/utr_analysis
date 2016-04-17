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
from utr.io import parse_input, create_build_dirs, combine_gff_results, output_coordinates
from utr.mapping import map_reads, filter_mapped_reads
from utr.util import setup_loggers
from utr.sequence import find_sequence, compute_coordinates

# Hide Biopython deprecation warnings
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

#--------------------------------------
# Global variables
#--------------------------------------

# Parse command-line arguments
args = parse_input()

# Get a list of sample ids

# For simplicity, the sample ID for a given file is assumed to be everything
# in the filename up to the mate pair specification.
input_regex = re.compile(r'^(.*/)?(?P<SAMPLE_ID>.*)_(R?[1-2])')

sample_ids = []

for filename in glob.glob(args.input_reads):
    sample_id = re.match(input_regex, filename).group('SAMPLE_ID')
    if sample_id not in sample_ids:
        sample_ids.append(sample_id)

# Create build directories
build_dirs = create_build_dirs(args, sample_ids)

# Setup global and task-specific loggers
loggers = setup_loggers(args.build_directory, build_dirs, sample_ids)

#-----------------------------------------------------------------------------
# RUFFUS TASKS
#
# Overview
# --------
# Below are the ruffus tasks associated with the UTR analysis pipeline. The
# first few tasks (0-2) are "common" tasks that are done regardless of the
# feature of interest (SL or Poly(A)). For example -- mapping against a
# "host" (e.g. host) genome to remove unrelated reads, and removing reads
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
    """check for bowtie2 indices and create if needed"""
    # check index for target species
    BASENAME_IDX = 0
    genome1 = os.path.splitext(args.target_genome)[BASENAME_IDX]

    # if index exists, stop here
    if os.path.exists('%s.1.bt2' % genome1):
        return

    # otherwise, create bowtie2 index
    logging.info("# Building bowtie2 index for %s" % args.target_genome)
    bowtie_cmd = (['bowtie2-build', args.target_genome, genome1])
    logging.info("# Command:\n" + " ".join(bowtie_cmd))
    ret = subprocess.call(bowtie_cmd)

    # stop here if we are only mapping to one genome
    if not args.host_genome:
        return

    # check index for host species
    genome2 = os.path.splitext(args.host_genome)[BASENAME_IDX]

    # if index exists, stop here
    if os.path.exists('%s.1.bt2' % genome2):
        return

    # otherwise, create bowtie2 index
    logging.info("# Building bowtie2 index for %s" % args.host_genome)
    bowtie_cmd = (['bowtie2-build', args.host_genome, genome2])
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
    BASENAME_IDX = 0
    target_genome = os.path.splitext(args.target_genome)[BASENAME_IDX]

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
    if not args.host_genome:
        return

    # Next, check filter genome
    host_genome = os.path.splitext(args.host_genome)[BASENAME_IDX]

    # if index exists, continue to next genome
    if os.path.exists('%s.fa' % host_genome):
        pass
    elif os.path.exists('%s.fasta' % host_genome):
        # if .fasta file exists, but not .fa, create a symlink
        logging.info("# Creating symlink to %s for Tophat" % args.host_genome)
        os.symlink(host_genome + '.fasta', host_genome + '.fa')
    else:
        raise IOError("Missing filter genome file")

#-----------------------------------------------------------------------------
# Step 1: Filter host reads (Optional)
#
# Before looking for spliced leader (SL) and Poly(A) reads, we will first
# remove any reads which come from an unrelated species such as host cells used
# to culture the target species.
#-----------------------------------------------------------------------------
@follows(check_for_bowtie_indices)
@follows(check_genome_fastas)
@collate(args.input_reads,
         formatter(r'^(.*/)?(?P<SAMPLE_ID>.*)_(R?[1-2])(.*)?\.fastq(?P<EXT>\.gz)?'),
         ['%s/{SAMPLE_ID[0]}/fastq/{SAMPLE_ID[0]}_host_reads_removed.1.fastq{EXT[0]}' % build_dirs['shared'],
          '%s/{SAMPLE_ID[0]}/fastq/{SAMPLE_ID[0]}_host_reads_removed.2.fastq{EXT[0]}' % build_dirs['shared']],
          '{SAMPLE_ID[0]}',
          '{EXT[0]}')
def filter_host_reads(input_reads, output_reads, sample_id, ext):
    # If we are only mapping to a single genome, we can just create
    # symlinks to the original input files and skip this step
    if not args.host_genome:
        # Determine which input corresponds to R1
        #if '1.fastq' in input_reads[0]:
        #    r1 = input_reads[0]
        #else:
        #    r1 = input_reads[1]

        #r2 = r1.replace('1.fastq', '2.fastq')

        ## Create symlinks
        #os.symlink(r1, output_file)
        #os.symlink(r2, output_file.replace('1.fastq', '2.fastq'))
        for read_num in [0,1]:
            if not os.path.exists(output_reads[read_num]):
                os.symlink(input_reads[read_num], output_reads[read_num])
        logging.info("# Skipping host read filtering step.")
        return

    logging.info("# Removing host reads.")

    # output directories
    tophat_dir = os.path.join(build_dirs['shared'], sample_id,
                             'tophat', 'mapped_to_host')
    fastq_dir = os.path.join(build_dirs['shared'], sample_id, 'fastq')

    # output fastq filepaths
    output_fastq = os.path.join(
        fastq_dir, "%s_host_reads_removed.fastq%s" % (sample_id, ext))

    # check for gff
    gff = args.host_gff if args.host_gff else None

    # map reads and remove hits
    filter_mapped_reads(input_reads[0], input_reads[1], args.host_genome,
                        tophat_dir, output_fastq, logging, gff=gff,
                        num_threads_tophat=args.num_threads_tophat)
    logging.info("# Finished removing host reads.")

#-----------------------------------------------------------------------------
# Step 2: Filter genomic reads
#
# The next step is to remove reads which map to the target genome before any
# trimming. These reads likely represent actual features in the genome and not
# the capped or polyadenylated reads we are interested in.
#-----------------------------------------------------------------------------
@subdivide(filter_host_reads,
           formatter(r'^(.*/)?(?P<SAMPLE_ID>.*)_host_reads_removed.(?P<READ_NUM>R?[1-2])\.fastq(?P<EXT>\.gz)?'),
           r'%s/{SAMPLE_ID[0]}/fastq/{SAMPLE_ID[0]}_genomic_reads_removed.*.fastq{EXT[0]}',
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}',
           '{EXT[0]}')
def filter_genomic_reads(input_files, output_files, sample_id, read_num, ext):
    # Clean-up any files from previous runs
    # for x in output_files:
        # logging.info("# Cleaning up pre-existing file: %s" % x)
        # if os.path.exists(x):
            # os.unlink(x)
    logging.info("# Removing genomic reads.")

    # output tophat and fastq directories
    tophat_dir = os.path.join(build_dirs['shared'], sample_id,
                              'tophat', 'mapped_to_target_untrimmed')
    output_fastq_dir = os.path.join(build_dirs['shared'], sample_id, 'fastq')

    # output fastq filepaths
    output_fastq = os.path.join(
        output_fastq_dir, "%s_genomic_reads_removed.fastq%s" % (sample_id, ext))

    # map reads and remove hits
    # Initially we will keep all (unmapped) reads which differ from genome by
    # at least one base. Later on we can be more restictive in our filtering
    # to make sure we aren't getting spurious hits.
    filter_mapped_reads(input_files[0], input_files[1],
                        args.target_genome, tophat_dir, output_fastq,
                        logging, read_mismatches=1, gff=args.target_gff,
                        num_threads_tophat=args.num_threads_tophat)
    logging.info("# Finished removing genomic reads.")

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
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)_genomic_reads_removed.(?P<READ_NUM>R?[12]).fastq(?P<EXT>\.gz)?'),
           '%s/{SAMPLE_ID[0]}/fastq/unfiltered/{SAMPLE_ID[0]}_{READ_NUM[0]}_1_sl_trimmed.fastq{EXT[0]}' % build_dirs['sl'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def find_sl_reads(input_reads, output_file, sample_id, read_num):
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

    log_handle = loggers[sample_id]['sl'][read_num]

    find_sequence(input_reads, 'sl', sl_filter, sl_regex, build_dirs['sl'],
                  sample_id, read_num, 
                  args.minimum_trimmed_length, args.max_dist_from_edge, 
                  log_handle)

#-----------------------------------------------------------------------------
# Step 4: Map trimmed reads
#
# Trims the matched sequences from reads and map to genome. For reads where
# the matched sequence comes from a trans-splicing or polyadenylation event,
# the location of the mapped trimmed read is where the addition took place.
#-----------------------------------------------------------------------------
@transform(find_sl_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)_(?P<READ_NUM>R?[12])_1_sl_trimmed.fastq(?P<EXT>\.gz)?'),
           '%s/{SAMPLE_ID[0]}/tophat/{READ_NUM[0]}_filtered_trimmed/accepted_hits.bam' % build_dirs['sl'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def map_sl_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered spliced-leader containing reads back to the genome"""
    log_handle = loggers[sample_id]['sl'][read_num]
    map_reads('sl', build_dirs['sl'], sample_id, read_num, 
              args.num_threads_tophat, log_handle)

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
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)/tophat/(?P<READ_NUM>R?[12])_filtered_trimmed/accepted_hits.bam'),
           '%s/{SAMPLE_ID[0]}/results/sl_coordinates_{READ_NUM[0]}.gff' % build_dirs['sl'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def compute_sl_coordinates(input_file, output_file, sample_id, read_num):
    logging.info("# Computing coordinates for mapped trans-splicing events "
                 "(%s)" % sample_id)
    log_handle = loggers[sample_id]['sl'][read_num]
    compute_coordinates(args.target_genome, 'sl', build_dirs['sl'], sample_id,
                        read_num, args.min_sl_length, 
                        args.minimum_differences, args.window_size, log_handle)
    logging.info("# Finished!")    

#-----------------------------------------------------------------------------
# Reverse SL analysis
#-----------------------------------------------------------------------------

#
# Reverse SL Step 1
#
@transform(filter_genomic_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)_genomic_reads_removed.(?P<READ_NUM>R?[12]).fastq(?P<EXT>\.gz)?'),
           '%s/{SAMPLE_ID[0]}/fastq/unfiltered/{SAMPLE_ID[0]}_{READ_NUM[0]}_1_rsl_trimmed.fastq{EXT[0]}' % build_dirs['rsl'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def find_rsl_reads(input_reads, output_file, sample_id, read_num):
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

    logging.info("# find_rsl_reads (%s %s)" % (sample_id, read_num))

    log_handle = loggers[sample_id]['rsl'][read_num]

    find_sequence(input_reads, 'rsl', rsl_filter, rsl_regex,
                  build_dirs['rsl'], sample_id, read_num, 
                  args.minimum_trimmed_length, args.max_dist_from_edge, 
                  log_handle,
                  trimmed_side='right', reverse=True)

#
# RSL Step 2
#
@transform(find_rsl_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)_(?P<READ_NUM>R?[12])_1_rsl_trimmed.fastq(?P<EXT>\.gz)?'),
           '%s/{SAMPLE_ID[0]}/tophat/{READ_NUM[0]}_filtered_trimmed/accepted_hits.bam' % build_dirs['rsl'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def map_rsl_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered poly-adenylated reads back to the genome"""
    log_handle = loggers[sample_id]['rsl'][read_num]
    map_reads('rsl', build_dirs['rsl'], sample_id, read_num,
              args.num_threads_tophat, log_handle)

#
# RSL Step 3
#
@transform(map_rsl_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)/tophat/(?P<READ_NUM>R?[12])_filtered_trimmed/accepted_hits.bam'),
           '%s/{SAMPLE_ID[0]}/results/rsl_coordinates_{READ_NUM[0]}.gff' % build_dirs['rsl'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def compute_rsl_coordinates(input_file, output_file, sample_id, read_num):
    logging.info("# Computing coordinates for mapped trans-splicing events "
                 "(%s, reverse)" % sample_id)
    log_handle = loggers[sample_id]['rsl'][read_num]
    compute_coordinates(args.target_genome, 'rsl', build_dirs['rsl'],
                        sample_id, read_num, args.min_sl_length, 
                        args.minimum_differences, args.window_size, log_handle)
    logging.info("# Finished!")    

#-----------------------------------------------------------------------------
# Poly(A) Analysis
#-----------------------------------------------------------------------------

#
# Poly(A) Step 1
#
@transform(filter_genomic_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)_genomic_reads_removed.(?P<READ_NUM>R?[12]).fastq(?P<EXT>\.gz)?'),
           '%s/{SAMPLE_ID[0]}/fastq/unfiltered/{SAMPLE_ID[0]}_{READ_NUM[0]}_1_polya_trimmed.fastq{EXT[0]}' % build_dirs['polya'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def find_polya_reads(input_reads, output_file, sample_id, read_num):
    """Matches reads with possible Poly(A) tail fragment"""
    polya_filter = 'A' * args.min_polya_length

    polya_regex = 'A{%d,}' % (args.min_polya_length)
    if args.exclude_internal_polya_matches:
        polya_regex = polya_regex + "$"

    logging.info("# find_polya_reads (%s %s)" % (sample_id, read_num))

    log_handle = loggers[sample_id]['polya'][read_num]

    find_sequence(input_reads, 'polya', polya_filter, polya_regex,
                  build_dirs['polya'], sample_id, read_num, 
                  args.minimum_trimmed_length, args.max_dist_from_edge, 
                  log_handle,
                  trimmed_side='right')

#
# Poly(A) Step 2
#
@transform(find_polya_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)_(?P<READ_NUM>R?[12])_1_polya_trimmed.fastq(?P<EXT>\.gz)?'),
           '%s/{SAMPLE_ID[0]}/tophat/{READ_NUM[0]}_filtered_trimmed/accepted_hits.bam' % build_dirs['polya'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def map_polya_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered poly-adenylated reads back to the genome"""
    log_handle = loggers[sample_id]['polya'][read_num]
    map_reads('polya', build_dirs['polya'], sample_id, read_num,
              args.num_threads_tophat, log_handle)

#
# Poly(A) Step 3
#
@transform(map_polya_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)/tophat/(?P<READ_NUM>R?[12])_filtered_trimmed/accepted_hits.bam'),
           '%s/{SAMPLE_ID[0]}/results/polya_coordinates_{READ_NUM[0]}.gff' % build_dirs['polya'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def compute_polya_coordinates(input_file, output_file, sample_id, read_num):
    logging.info("# Computing coordinates for mapped polyadenylation sites "
                 "(%s)" % sample_id)
    log_handle = loggers[sample_id]['polya'][read_num]
    compute_coordinates(args.target_genome, 'polya', build_dirs['polya'],
                        sample_id, read_num, args.min_polya_length,
                        args.minimum_differences, args.window_size, log_handle)
    logging.info("# Finished!")    

#-----------------------------------------------------------------------------
# Poly(T) Analysis
#-----------------------------------------------------------------------------

#
# Poly(T) Step 1
#
@transform(filter_genomic_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)_genomic_reads_removed.(?P<READ_NUM>R?[12]).fastq(?P<EXT>\.gz)?'),
           '%s/{SAMPLE_ID[0]}/fastq/unfiltered/{SAMPLE_ID[0]}_{READ_NUM[0]}_1_polyt_trimmed.fastq{EXT[0]}' % build_dirs['polyt'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def find_polyt_reads(input_reads, output_file, sample_id, read_num):
    """Matches reads with possible Poly(T) tail fragment"""
    # Match reads with at least n T's at the beginning of the read; For now
    # we will always require matches to be at the beginning of the read.
    polyt_filter = 'T' * args.min_polya_length

    polyt_regex = 'T{%d,}' % (args.min_polya_length)
    if args.exclude_internal_polya_matches:
        polyt_regex = '^' + polyt_regex

    logging.info("# find_polyt_reads (%s %s)" % (sample_id, read_num))

    log_handle = loggers[sample_id]['polyt'][read_num]
    
    find_sequence(input_reads, 'polyt', polyt_filter, polyt_regex,
                  build_dirs['polyt'], sample_id, read_num, 
                  args.minimum_trimmed_length, args.max_dist_from_edge, 
                  log_handle,
                  reverse=True)

#
# Poly(T) Step 2
#
@transform(find_polyt_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)_(?P<READ_NUM>R?[12])_1_polyt_trimmed.fastq(?P<EXT>\.gz)?'),
           '%s/{SAMPLE_ID[0]}/tophat/{READ_NUM[0]}_filtered_trimmed/accepted_hits.bam' % build_dirs['polyt'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def map_polyt_reads(input_file, output_file, sample_id, read_num):
    """Maps the filtered Poly(T) reads back to the genome"""
    log_handle = loggers[sample_id]['polyt'][read_num]
    map_reads('polyt', build_dirs['polyt'], sample_id, read_num,
              args.num_threads_tophat, log_handle)

#
# Poly(T) Step 3
#
@transform(map_polyt_reads,
           formatter(r'^(.*)/(?P<SAMPLE_ID>.+)/tophat/(?P<READ_NUM>R?[12])_filtered_trimmed/accepted_hits.bam'),
           '%s/{SAMPLE_ID[0]}/results/polyt_coordinates_{READ_NUM[0]}.gff' % build_dirs['polyt'],
           '{SAMPLE_ID[0]}',
           '{READ_NUM[0]}')
def compute_polyt_coordinates(input_file, output_file, sample_id, read_num):
    logging.info("# Computing coordinates for mapped polyadenylation sites "
                 "(%s, reverse)" % sample_id)
    log_handle = loggers[sample_id]['polyt'][read_num]
    compute_coordinates(args.target_genome, 'polyt', build_dirs['polyt'],
                        sample_id, read_num, args.min_polya_length,
                        args.minimum_differences, args.window_size, log_handle)
    logging.info("# Finished!")    

#-----------------------------------------------------------------------------
# Step 6: Combine coordinate output for multiple samples
#
# Next, we create a single GFF file using the information contained in the GFF
# files that were constructed for each sample.
#
# Coordinates with low coverage may also be filtered out at this step.
#-----------------------------------------------------------------------------
@merge([compute_sl_coordinates,
        compute_rsl_coordinates],
        # formatter(r'^(.*/)/(?P<SAMPLE_ID>.*)/results/(?P<FEATURE_TYPE>[a-z]+)_coordinates_(?P<READ_NUM>R?[1-2]).gff'),
        os.path.join(build_dirs['combined'], 'spliced_leader.gff'))
def combine_sl_results(input_files, output_file):
    # Combine spliced leader output
    logging.info("# Combining spliced leader coordinates output")

    sl_combined = combine_gff_results(input_files)

    # Save summary GFFs
    sl_outfile = os.path.join(build_dirs['combined'], 'spliced_leader.gff')
    output_coordinates(sl_combined, 'sl', sl_outfile, track_color='83,166,156')

@merge([compute_polya_coordinates,
        compute_polyt_coordinates],
        # formatter(r'^(.*/)/(?P<SAMPLE_ID>.*)/results/(?P<FEATURE_TYPE>[a-z]+)_coordinates_(?P<READ_NUM>R?[1-2]).gff'),
        os.path.join(build_dirs['combined'], 'polya.gff'))
def combine_polya_results(input_files, output_file):
    # Combine spliced leader output
    logging.info("# Combining spliced leader coordinates output")

    polya_combined = combine_gff_results(input_files)

    # Save summary GFFs
    polya_outfile = os.path.join(build_dirs['combined'], 'polya.gff')
    output_coordinates(polya_combined, 'polya', polya_outfile,
                       track_color='166,83,93')

#-----------------------------------------------------------------------------
# Run pipeline
#-----------------------------------------------------------------------------
if __name__ == "__main__":
    pipeline_run([combine_sl_results, combine_polya_results],
                 logger=logging.getLogger(''),
                 multiprocess=args.num_threads)
    pipeline_printout_graph("utr_analysis_flowchart.png", "png",
                            [combine_sl_results, combine_polya_results],
                            pipeline_name='Trypanosome UTR Analysis Pipeline')

