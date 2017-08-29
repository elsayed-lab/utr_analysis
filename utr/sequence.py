"""
Sequence and gene structure related functions:

    - find_sequence
    - compute_coordinates
    - is_inside_cds
    - find_closest_gene

"""
import csv
import distance
import gzip
import os
import pysam
import re
import cStringIO
from Bio import SeqIO
from io import output_coordinates, load_annotations, readfq
from util import num_lines, compress_str

def find_sequence(input_file, feature_name, sequence_filter, feature_regex,
                  build_dir, sample_id, read_num, minimum_trimmed_length,
                  max_dist_from_edge, log_handle):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing a specified sequence of interest.

    Arguments
    ---------
    input_file: str
        Filepath to a FASTQ file containing reads to scan.
    feature_name: str
        Type of feature being searched for; used in naming filing and
        directories and in choosing logs to write to. [sl|polya]
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
        Which of the mated reads should be scanned. [1|2]
    minimum_trimmed_length: int
        Minimum length of read allowed after matching feature is trimmed.
    max_dist_from_edge: int
        Maximum distance SL/Poly(A) feature can be from the edge of read.
    log_handle: logging.Handle
        Handler to use for logging.

    Output files
    ------------
    There are three possible sets of output files for this function depending
    on whether the input read comes from a mated-pair of reads or a single
    read, and whether (for the case of mated-pair reads) it is the left read
    or right read:

    1. Sequence found in R1
        *_1_1_xxx_untrimmed.fastq
        *_1_1_xxx_trimmed.fastq
        *_1_2.fastq
    2. Sequence found in R2
        *_2_2_xxx_untrimmed.fastq
        *_2_2_xxx_trimmed.fastq
        *_2_1.fastq
    """
    #--------------------------------------
    # FASTQ row indices
    #--------------------------------------
    ID_IDX = 0
    SEQUENCE_IDX = 1
    QUALITY_IDX = 2

    log_handle.info("# Processing %s" % os.path.basename(input_file))

    # determine whether regular expression in anchored
    anchored = (feature_regex.startswith("^") or feature_regex.endswith("$"))

    # list to keep track of potential matches
    matches = []

    # output filepaths
    output_base = '%s/%s/fastq/%s_%s_%s' % (
        build_dir, sample_id, sample_id, read_num, read_num[-1]
    )

    # determine compression type to use
    file_ext = os.path.splitext(input_file)[-1]

    # for uncompressed fastq files, we don't need the final extension
    if file_ext not in ['.gz', '.xz']:
        file_ext = ''

    output_untrimmed = "%s_%s_untrimmed.fastq%s" % (output_base, feature_name,
                                                     file_ext)
    output_trimmed = "%s_%s_trimmed.fastq%s" % (output_base, feature_name,
                                                file_ext)

    # Also keep track of match lengths which will be used for more rigorous
    # filtering when comparing to the genome sequence near where the read is
    # mapped.
    match_lengths_dir = os.path.join(build_dir, sample_id, 'results')
    output_lengths = "%s/match_lengths_%s.csv.gz" % (match_lengths_dir, read_num)
    match_lengths_fp = gzip.open(output_lengths, 'wb')

    # mated reads
    read_num_other = "1" if read_num == "2" else "2"
    input_file_mated = input_file.replace("." + read_num[-1],
                                          "." + read_num_other[-1])
    output_mated_reads = "%s_%s.fastq%s" % (output_base[:-2], read_num_other,
                                             file_ext)

    # compile regex
    read_regex = re.compile(feature_regex)

    # total number of reads
    num_reads = num_lines(input_file) / 4

    # Start sample log
    log_handle.info("# Scanning %d reads for %s" % (num_reads, feature_name))
    log_handle.info("# Using Regex pattern:\n %s" % feature_regex)

    # open output string buffer (will write to compressed file later)
    reads_trimmed = cStringIO.StringIO()
    reads_untrimmed = cStringIO.StringIO()
    mated_reads_buffer = cStringIO.StringIO()

    # Keep track of matched read IDs
    read_ids = []

    # Keep track of ways in which reads are filtered out
    num_filtered_no_seq_match = 0
    num_filtered_too_small = 0
    num_filtered_far_from_edge = 0

    # Find all reads containing the sequence of interest
    if file_ext == '.gz':
        fastq = gzip.open(input_file, 'rb')
        fastq_mated = gzip.open(input_file_mated, 'rb')
    elif file_ext == '.xz':
        import backports.lzma as lzma
        fastq = lzma.open(input_file, 'rb')
        fastq_mated = lzma.open(input_file_mated, 'rb')
    else:
        fastq = open(input_file, 'r')
        fastq_mated = open(input_file_mated, 'r')

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
            num_filtered_no_seq_match += 1
            continue

        # check for match

        # When looking for internal matches, there may be multiple hits. Choose
        # the one that is closest to the edge of the read where the feature is
        # expected to be found.
        try:
            # for polya, reverse sequence to to find match closest to right
            # side
            if (feature_name == 'polya' and (not anchored)):
                match = re.search(read_regex, read[SEQUENCE_IDX][::-1])
                match_start = len(read[SEQUENCE_IDX]) - match.end() 
                match_end = len(read[SEQUENCE_IDX]) - match.start() 
            else:
                match = re.search(read_regex, read[SEQUENCE_IDX])

                # for anchored reads, its possible that a read passes the
                # quick filter check but the regex does not match
                if match is None:
                    num_filtered_no_seq_match += 1
                    continue

                match_start = match.start()
                match_end = match.end()
        except:
            import pdb; pdb.set_trace();

        # match length
        match_length = match.end() - match.start()

        # For SL sequence, trim everything up to end of match
        if feature_name == 'sl':
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

        # skip reads that are less than the required amount after trimming
        if len(trimmed_read[SEQUENCE_IDX]) < minimum_trimmed_length:
            num_filtered_too_small += 1
            continue

        # length of portion trimmed off
        trimmed_part_length = (len(read[SEQUENCE_IDX]) - 
                               len(trimmed_read[SEQUENCE_IDX]))

        # for internal matches, skip reads where match is not close enough to
        # the edge of the read
        if (trimmed_part_length - match_length) > max_dist_from_edge:
            num_filtered_far_from_edge += 1
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
    log_handle.info("# Found %d reads with possible %s fragment" %
             (len(read_ids), feature_name))
    log_handle.info(
        "# Excluded %d reads with no feature matches."
        % num_filtered_no_seq_match)
    log_handle.info(
        "# Excluded %d reads which were too short after trimming."
        % num_filtered_too_small)
    log_handle.info(
        "# Excluded %d reads with matched feature too far from read edge."
        % num_filtered_far_from_edge)

    # Create output directory
    output_dir = os.path.dirname(output_base)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, mode=0o755)

    # write trimmed and untrimmed reads to fastq.gz / fastq.xz
    compress_str(output_trimmed, reads_trimmed)
    compress_str(output_untrimmed, reads_untrimmed)
    compress_str(output_mated_reads, mated_reads_buffer)

    # clean up
    fastq.close()
    fastq_mated.close()
    reads_trimmed.close()
    reads_untrimmed.close()
    mated_reads_buffer.close()
    match_lengths_fp.close()

    log_handle.info("# Finished processing %s" % os.path.basename(input_file))

def compute_coordinates(target_genome, target_gff, feature_name, build_dir, 
                        shared_build_dir, sample_id,
                        read_num, min_feature_length, minimum_differences, 
                        window_size, log_handle):
    """
    Outputs coordinates and frequencies of putative UTR features to a GFF file.

    Sequences manipulated below include:

        untrimmed_seq       untrimmed read sequence
        trimmed_seq         sequence of read after trimming
        trimmed_part        sequence of the portion of the read trimmed off
        trimmed_genome_seq  sequence at the genome location corresponding to
                            trimmed portion of the read
        matched_seq         SL/Poly(A) sequence matched in previous step
        matched_genome_seq  sequence at the genome location corresponding to
                            the matched portion of the read

    Note that, unless unanchored matching is used and the match is
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
    genome = SeqIO.parse(target_genome, 'fasta')
    chr_sequences = {x.id:x for x in genome}

    # Get chromosomes/contigs from GFF file
    chromosomes = load_annotations(target_gff)
    
    # Create a dictionary to keep track of the coordinates.
    results = {}

    # keep track of how many reads were found in the expected location
    num_good = 0

    # Keep track of ways in which reads are filtered out
    num_filtered_too_small = 0
    num_filtered_near_chr_edge = 0
    num_filtered_matches_genome_seq = 0
    num_filtered_similar_to_genome_seq = 0
    num_filtered_no_nearby_genes = 0
    num_filtered_inside_cds = 0

    # DEBUGGING: checking for one-off error when determining genome seq
    neg3 = 0
    neg2 = 0
    neg1 = 0
    zero = 0
    pos1 = 0
    pos2 = 0
    pos3 = 0
    pos_strand = 0
    neg_strand = 0

    pos_strand_late = 0
    neg_strand_late = 0
    pos_strand_late2 = 0
    neg_strand_late2 = 0
    pos_strand_late3 = 0
    neg_strand_late3 = 0
    pos_strand_late4 = 0
    neg_strand_late4 = 0

    # Create output CSV writer for hits that are not near any genes
    no_nearby_genes = csv.writer(
        open('%s/no_nearby_genes_%s.csv' % (output_dir, read_num), 'w')
    )
    no_nearby_genes.writerow(['read_id', 'chromosome', 'strand', 'position'])

    # Create output CSV writer for hits that are found inside known CDS's
    inside_cds = csv.writer(
        open('%s/inside_cds_%s.csv' % (output_dir, read_num), 'w')
    )
    inside_cds.writerow(['gene_id', 'chromosome', 'strand', 'position'])

    # In addition to saving the coordinates as a GFF file, we will also write a
    # CSV file which contains entries for each read used. This can be useful
    # for debugging/tracking down the origin of a particular coordinate.
    sample_csv_writer = csv.writer(
        open('%s/matched_reads_%s.csv' % (output_dir, read_num), 'w')
    )
    sample_csv_writer.writerow(['read_id', 'gene_id', 'chromosome',
        'read_strand', 'acceptor_site_side', 'debug_match_loc',
        'trimmed_start', 'trimmed_stop',
        'untrimmed_seq', 'trimmed_part', 'trimmed_genome_seq', 'matched_seq',
        'matched_genome_seq', 'acceptor_site'])

    # Load the untrimmed reads in order to to determine original read lengths
    input_bam_untrimmed = '%s/%s/tophat/mapped_to_target_untrimmed/unmapped_sorted.bam' % (
        shared_build_dir, sample_id
    )
    sam_untrimmed = pysam.Samfile(input_bam_untrimmed, 'rb')

    untrimmed_reads = {}

    for read in sam_untrimmed:
        rnum = '1' if read.is_read1 else '2'
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

    # Keep track of reads skipped due to coming from the wrong side of PE
    num_skipped = 0

    # Keep track of matched read IDs
    read_ids = []

    # Get coordinate and strand for each read in bam file
    for i, read in enumerate(sam):
        # Get read where feature sequence was found
        if ((read.is_read1 and read_num == '2') or
            (read.is_read2 and read_num == '1')):
            num_skipped += 1
            continue

        read_ids.append(read.qname)

        # Chromosome and strand where the read was mapped
        chromosome = sam.getrname(read.tid)
        strand = "-" if read.is_reverse else "+"

        if strand == '+':
            pos_strand = pos_strand + 1
        else:
            neg_strand = neg_strand + 1

        # temp
        if i % 1000 == 0:
            print("Processing read %d" % i)

        # Length of matched SL suffix or number of A's/T's
        untrimmed_read = untrimmed_reads[read.qname][read_num]

        # get match length (actual portion matched)
        match_length = match_lengths[read.qname]

        # if trimmed read mapped to reverse strand, then polya/acceptor site
        # will be on the opposite side of where it was originally located
        if ((strand == '+' and feature_name == 'sl') or 
            (strand == '-' and feature_name == 'polya')):
            acceptor_site_side = 'left'
        else:
            acceptor_site_side = 'right'

        # Acceptor site position
        # Note: csamtools.AlignedRead positions are 0-indexed
        if acceptor_site_side == 'left':
            acceptor_site = read.pos
        else:
            acceptor_site = read.pos + read.rlen + 1

        # If read is mapped too close to end of chromosome, skip it
        # left end
        if read.pos < (min_feature_length - 1):
            num_filtered_near_chr_edge += 1
            continue
        # right end
        elif ((read.pos + read.rlen) > (len(chromosomes[chromosome]) - min_feature_length)):
            num_filtered_near_chr_edge += 1
            continue

        # determine the length of the trimmed portion of read
        trimmed_part_length = get_trimmed_part_len(acceptor_site_side,
                                                   untrimmed_read.rlen,
                                                   read.rlen, read.pos, 
                                                   len(chr_sequences[chromosome]))

        # get the genomic sequences corresponding to the trimmed portion of
        # the read and the matched portion of the read
        trimmed_genome_seq = get_trimmed_genome_seq(acceptor_site_side,
                                                    read.pos, read.rlen, 
                                                    trimmed_part_length,
                                                    chr_sequences[chromosome],
                                                    strand)

        # read orientation has been normalized so we can check the same
        # location and direction for both strands
        matched_genome_seq = str(trimmed_genome_seq[-match_length:])

        # Get sequence of the trimmed portion of the read and the sequence that
        # matched feature of interest
        if feature_name == 'sl':
            # get trimmed portion from left of read
            trimmed_part = str(untrimmed_read.seq[:trimmed_part_length])
            matched_seq = trimmed_part[-match_length:]
        else:
            trimmed_part = str(untrimmed_read.seq[-trimmed_part_length:])
            matched_seq = trimmed_part[:match_length]

        # TEMP DEBUGGING
        if acceptor_site_side == 'left': 
            start = read.pos - trimmed_part_length
            end = read.pos
        else:
            start = read.pos + read.rlen
            end = read.pos + read.rlen + trimmed_part_length

        debug_match_loc = 'none'

        if feature_name == 'sl':
            if strand == '+':
                if str(chr_sequences[chromosome][start:end].seq[:match_length]).endswith('TTG'):
                    zero += 1
                    debug_match_loc = 'expected'
                elif str(chr_sequences[chromosome][start-3:end-3].seq[:match_length]).endswith('TTG'):
                    neg3 += 1
                    debug_match_loc = 'neg3'
                elif str(chr_sequences[chromosome][start-2:end-2].seq[:match_length]).endswith('TTG'):
                    neg2 += 1
                    debug_match_loc = 'neg2'
                elif str(chr_sequences[chromosome][start-1:end-1].seq[:match_length]).endswith('TTG'):
                    neg1 += 1
                    debug_match_loc = 'neg1'
                elif str(chr_sequences[chromosome][start+1:end+1].seq[:match_length]).endswith('TTG'):
                    debug_match_loc = 'pos1'
                    pos1 += 1
                elif str(chr_sequences[chromosome][start+2:end+2].seq[:match_length]).endswith('TTG'):
                    debug_match_loc = 'pos2'
                    pos2 += 1
                elif str(chr_sequences[chromosome][start+3:end+3].seq[:match_length]).endswith('TTG'):
                    pos3 += 1
                    debug_match_loc = 'pos3'
            # sl, neg
            else:
                if str(chr_sequences[chromosome][start:end].seq[:match_length].reverse_complement()).endswith('TTG'):
                    debug_match_loc = 'expected'
                    zero += 1
                elif str(chr_sequences[chromosome][start-3:end-3].seq[:match_length].reverse_complement()).endswith('TTG'):
                    debug_match_loc = 'neg3'
                    neg3 += 1
                elif str(chr_sequences[chromosome][start-2:end-2].seq[:match_length].reverse_complement()).endswith('TTG'):
                    debug_match_loc = 'neg2'
                    neg2 += 1
                elif str(chr_sequences[chromosome][start-1:end-1].seq[:match_length].reverse_complement()).endswith('TTG'):
                    debug_match_loc = 'neg1'
                    neg1 += 1
                elif str(chr_sequences[chromosome][start+1:end+1].seq[:match_length].reverse_complement()).endswith('TTG'):
                    debug_match_loc = 'pos1'
                    pos1 += 1
                elif str(chr_sequences[chromosome][start+2:end+2].seq[:match_length].reverse_complement()).endswith('TTG'):
                    debug_match_loc = 'pos2'
                    pos2 += 1
                elif str(chr_sequences[chromosome][start+3:end+3].seq[:match_length].reverse_complement()).endswith('TTG'):
                    debug_match_loc = 'pos3'
                    pos3 += 1
        else:
            if strand == '+':
                if str(chr_sequences[chromosome][start:end].seq[:match_length]).endswith('AAA'):
                    zero += 1
                    debug_match_loc = 'expected'
                elif str(chr_sequences[chromosome][start-3:end-3].seq[:match_length]).endswith('AAA'):
                    neg3 += 1
                    debug_match_loc = 'neg3'
                elif str(chr_sequences[chromosome][start-2:end-2].seq[:match_length]).endswith('AAA'):
                    neg2 += 1
                    debug_match_loc = 'neg2'
                elif str(chr_sequences[chromosome][start-1:end-1].seq[:match_length]).endswith('AAA'):
                    neg1 += 1
                    debug_match_loc = 'neg1'
                elif str(chr_sequences[chromosome][start+1:end+1].seq[:match_length]).endswith('AAA'):
                    debug_match_loc = 'pos1'
                    pos1 += 1
                elif str(chr_sequences[chromosome][start+2:end+2].seq[:match_length]).endswith('AAA'):
                    debug_match_loc = 'pos2'
                    pos2 += 1
                elif str(chr_sequences[chromosome][start+3:end+3].seq[:match_length]).endswith('AAA'):
                    pos3 += 1
                    debug_match_loc = 'pos3'
            else:
                if str(chr_sequences[chromosome][start:end].seq[:match_length]).endswith('TTT'):
                    zero += 1
                    debug_match_loc = 'expected'
                elif str(chr_sequences[chromosome][start-3:end-3].seq[:match_length]).endswith('TTT'):
                    neg3 += 1
                    debug_match_loc = 'neg3'
                elif str(chr_sequences[chromosome][start-2:end-2].seq[:match_length]).endswith('TTT'):
                    neg2 += 1
                    debug_match_loc = 'neg2'
                elif str(chr_sequences[chromosome][start-1:end-1].seq[:match_length]).endswith('TTT'):
                    neg1 += 1
                    debug_match_loc = 'neg1'
                elif str(chr_sequences[chromosome][start+1:end+1].seq[:match_length]).endswith('TTT'):
                    debug_match_loc = 'pos1'
                    pos1 += 1
                elif str(chr_sequences[chromosome][start+2:end+2].seq[:match_length]).endswith('TTT'):
                    debug_match_loc = 'pos2'
                    pos2 += 1
                elif str(chr_sequences[chromosome][start+3:end+3].seq[:match_length]).endswith('TTT'):
                    pos3 += 1
                    debug_match_loc = 'pos3'

        # Check to see that match differs from genome sequence (quick)
        if (matched_seq == matched_genome_seq):
            num_filtered_matches_genome_seq += 1
            continue

        # Note: for real acceptor sites, strand = gene_strand

        # For Poly(A) tail, extend detected read to include any A's present in
        # the genome that were trimmed off
        if feature_name == 'polya':
            match = re.search('^A*', matched_genome_seq)
            overlap_length = match.end() - match.start()

            if overlap_length > 0:
                matched_seq = matched_seq[overlap_length:]
                matched_genome_seq = matched_genome_seq[overlap_length:]

                if strand == '+':
                    # give back some A's and update the relevant sequences
                    #matched_genome_seq = matched_genome_seq[overlap_length:]
                    acceptor_site = acceptor_site + overlap_length
                elif strand == '-':
                    # give back some A's and update the relevant sequences
                    #matched_genome_seq = matched_genome_seq[-overlap_length:]
                    acceptor_site = acceptor_site - overlap_length

        # 
        # More filtering
        #
        # Check to see if feature is still sufficiently long;
        if len(matched_seq) < min_feature_length:
            num_filtered_too_small += 1
            continue

        # Check to see that match (still) differs from genome sequence
        if (matched_seq == matched_genome_seq):
            num_filtered_matches_genome_seq += 1
            continue

        if strand == '+':
            pos_strand_late += 1
        else:
            neg_strand_late += 1

        # Make sure that for the region matched, there are at least n
        # differences between the match and the corresponding genomic
        # sequence (slower)
        seq_dist = distance.hamming(matched_seq, matched_genome_seq)
        if seq_dist < minimum_differences:
            num_filtered_similar_to_genome_seq += 1
            continue

        if strand == '+':
            pos_strand_late2 += 1
        else:
            neg_strand_late2 += 1

        # Check to make sure the acceptor site does not fall within
        # a known CDS: if it does, save to a separate file to look at later
        cds_id = is_inside_cds(chromosomes[chromosome], acceptor_site, window_size)

        if cds_id != "":
            inside_cds.writerow([cds_id, chromosome, strand, acceptor_site])
            num_filtered_inside_cds += 1
            continue

        if strand == '+':
            pos_strand_late3 += 1
        else:
            neg_strand_late3 += 1

        # Find nearest gene
        gene = find_closest_gene(chromosomes[chromosome], feature_name,
                                 acceptor_site, acceptor_site_side, window_size)

        # If no nearby genes were found, stop here
        if gene is None:
            no_nearby_genes.writerow(
                [read.qname, chromosome, strand, acceptor_site])
            num_filtered_no_nearby_genes += 1
            continue

        if strand == '+':
            pos_strand_late4 += 1
        else:
            neg_strand_late4 += 1

        num_good = num_good + 1

        # Add to output dictionary
        if not chromosome in results:
            results[chromosome] = {}
        if not gene['id'] in results[chromosome]:
            results[chromosome][gene['id']] = {}

        # Add entry to sample output csv
        sample_csv_writer.writerow([
            read.qname, gene['id'], chromosome, strand, 
            acceptor_site_side, debug_match_loc,
            read.pos, read.pos + read.rlen, untrimmed_read.seq, trimmed_part,
            str(trimmed_genome_seq), matched_seq, matched_genome_seq, 
            acceptor_site
        ])

        #
        # Increment site count and save distance from gene
        #
        # count: number of reads supporting a given site
        # distance: distance between the site and the nearest CDS
        # seq_dist: hamming distance between matched sl/poly(A) sequence
        #           and the sequence at the genome for that location
        # match_len: length of SL/Poly(A) sequence matched
        # strand: strand of mapped site
        # description: gene description from original GFF
        #
        if not acceptor_site in results[chromosome][gene['id']]:
            results[chromosome][gene['id']][acceptor_site] = {
                "count": 1,
                "distance": gene['distance'],
                "seq_dist": seq_dist,
                "match_len": len(matched_seq),
                "strand": strand,
                "description": gene['description']
            }
        else:
            results[chromosome][gene['id']][acceptor_site]['count'] += 1
            results[chromosome][gene['id']][acceptor_site]['seq_dist'] += seq_dist
            results[chromosome][gene['id']][acceptor_site]['match_len'] += len(matched_seq)

    # total number of reads checked
    total_reads = float(i) - num_skipped

    # record number of good and bad reads
    log_handle.info("#########################################################")
    log_handle.info("#")
    log_handle.info("# Scanned a total of %d reads for %s fragments:" %
                    (i, feature_name))
    log_handle.info("#")
    log_handle.info("#    - Matches: %d (%0.2f%%)" % (num_good, (num_good / total_reads) * 100))
    log_handle.info("#    - Excluded (too short): %d (%0.2f%%)" % 
                    (num_filtered_too_small, (num_filtered_too_small / total_reads) * 100))
    log_handle.info(
        "#    - Excluded (near chromosome boundaries): %d (%0.2f%%)" % 
        (num_filtered_near_chr_edge, (num_filtered_near_chr_edge / total_reads) * 100))
    log_handle.info(
        "#    - Excluded (matches genome sequence): %d (%0.2f%%)" %
        (num_filtered_matches_genome_seq, (num_filtered_matches_genome_seq / total_reads) * 100))
    log_handle.info(
        "#    - Excluded (simmilar to genome sequence): %d (%0.2f%%)" % 
        (num_filtered_similar_to_genome_seq, (num_filtered_similar_to_genome_seq / total_reads) * 100))
    log_handle.info(
        "#    - Excluded (inside known CDS): %d (%0.2f%%)"  % 
        (num_filtered_inside_cds, (num_filtered_inside_cds / total_reads) * 100))
    log_handle.info(
        "#    - Excluded (to far from known CDS's): %d (%0.2f%%)" % 
        (num_filtered_no_nearby_genes, (num_filtered_no_nearby_genes / total_reads) * 100))
    log_handle.info("#")
    log_handle.info("# DEBUGGING:")
    log_handle.info("#")
    log_handle.info("#    - neg3: %d" % neg3)
    log_handle.info("#    - neg2: %d" % neg2)
    log_handle.info("#    - neg1: %d" % neg1)
    log_handle.info("#    - zero: %d" % zero)
    log_handle.info("#    - pos1: %d" % pos1)
    log_handle.info("#    - pos2: %d" % pos2)
    log_handle.info("#    - pos3: %d" % pos3)
    log_handle.info("#")
    log_handle.info("#    - + strand: %d (late: %d / %d / %d / %d)" % (pos_strand, pos_strand_late, pos_strand_late2, pos_strand_late3, pos_strand_late4))
    log_handle.info("#    - - strand: %d (late: %d / %d / %d / %d)" % (neg_strand, neg_strand_late, neg_strand_late2, neg_strand_late3, neg_strand_late4))
    log_handle.info("#")
    log_handle.info("#########################################################")

    # Output coordinates as a GFF file
    output_filepath = '%s/%s_coordinates_%s.gff' % (
        output_dir, feature_name, read_num
    )

    output_coordinates(results, feature_name, output_filepath, target_gff)

def get_trimmed_part_len(acceptor_site_side, untrimmed_read_len,
                         read_len, read_pos, chr_len):
    """
    Returns the length of the region of the read that was trimmed off.

    In most cases, this is simply the difference between the untrimmed read
    length and and trimmed read length.

    Near chromosome boundaries, however, the trimmed part length may need to be
    truncated.

    NOTE: this may be larger than the SL sequence length if unanchored
    matches are allowed
    """
    trimmed_part_length = untrimmed_read_len - read_len

    if acceptor_site_side == 'left': 
        # Shorten read if left end of read was mapped near end of
        # chromosome
        trimmed_part_length = min(min(trimmed_part_length, read_pos),
                                    chr_len - read_pos)
        
    # Genome sequence just downstream of mapped location
    else:
        # Shorten if right end of read was mapped near end of chromosome
        trimmed_part_length = min(min(trimmed_part_length, read_pos),
                                    chr_len - (read_pos + read_len))

    return trimmed_part_length


def get_trimmed_genome_seq(acceptor_site_side, read_pos, read_len,
                            trimmed_part_length, chr_seq, strand):
    """Determines the genomic sequence corresponding to the portion of
    the read which was trimmed off."""
    # Genome sequence just upstream of mapped location
    if acceptor_site_side == 'left': 
        # Grab region just before mapped read
        start = read_pos - trimmed_part_length
        end = read_pos
        # end = read_pos + 1

        trimmed_genome_seq = chr_seq[start:end].seq
    # Genome sequence just downstream of mapped location
    else:
        # Grab region just after mapped untrimmed read
        start = read_pos + read_len
        end = read_pos + read_len + trimmed_part_length

        # TESTING 2016/12/02
        # start = read_pos + read_len - 1
        # end = read_pos + read_len + trimmed_part_length - 1
        trimmed_genome_seq = chr_seq[start:end].seq

    # For reads mapped to negative strand, take complement
    if strand == '+':
        trimmed_genome_seq = trimmed_genome_seq
    else:
        trimmed_genome_seq = trimmed_genome_seq.reverse_complement()

    return trimmed_genome_seq

def is_inside_cds(chromosome, location, window_size):
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
    str:
        ID of CDS if found to be inside of one, otherwise ""

    References
    ----------
    http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html
    """
    # Number of bases before or after feature
    half_window = int(window_size / 2)

    # Extract region of sequence surrounding the putative feature of
    # interest
    nearby = chromosome[max(location - half_window, 0):location + half_window]

    # determine relative location of acceptor_site
    if location < half_window:
        # at left-end, relative site is the actual one
        relative_location = location
    else:
        # otherwise, window will be centered around the acceptor site, or, if
        # near the right-end of the chromosome, the acceptor site will fall
        # at the half-window mark
        #relative_location = len(nearby) / 2
        relative_location = half_window

    # scan all genes near putative feature
    for feature in nearby.features:
        if ((feature.location.start <= relative_location) and
            (feature.location.end >= relative_location)):
            return feature.id

    return ""

def find_closest_gene(chromosome, feature_name, location, acceptor_site_side,
                      window_size):
    """
    Finds the closest gene to a specified location that is in the expected
    orientation.
    """
    # chromosome boundary
    ch_end = len(chromosome)

    # 1. Get genes within +/- window_size/2 bases of location
    half_win = window_size / 2

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

        # For SL, look at gene start locations and for Poly(A) look
        # at where each gene ends.
        dist = abs(gene.location.start - feature_location)

        # For Poly(A) look at gene endings
        if dist < closest_dist:
            closest_index = i
            closest_gene = gene
            closest_dist = dist

    # No genes found in correct orientation
    if closest_gene is None:
        return None

    # Make sure strand of gene is appropriate for the feature and
    # orientation
    if acceptor_site_side == 'left':
        if ((closest_gene.strand == -1 and feature_name == 'sl') or 
            (closest_gene.strand ==  1 and feature_name == 'polya')):
            return None
    elif acceptor_site_side == 'right':
        if ((closest_gene.strand == -1 and feature_name == 'polya') or 
            (closest_gene.strand ==  1 and feature_name == 'sl')):
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
        'id': closest_gene.id,
        'description': gene_description,
        'distance': closest_dist
    }
