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
import StringIO
from Bio import SeqIO
from io import output_coordinates

def find_sequence(input_file, feature_name, sequence_filter, feature_regex,
                  build_dir, sample_id, read_num, minimum_trimmed_length,
                  max_dist_from_edge, log_handle,
                  trimmed_side='left', reverse=False):
    """
    Loads a collection of RNA-Seq reads and filters the reads so as to only
    return those containing a specified sequence of interest.

    Arguments
    ---------
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
    minimum_trimmed_length: int
        Minimum length of read allowed after matching feature is trimmed.
    max_dist_from_edge: int
        Maximum distance SL/Poly(A) feature can be from the edge of read.
    log_handle: logging.Handle
        Handler to use for logging.
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
    log_handle.info("# Scanning %d reads for %s" % (num_reads, feature_name))
    log_handle.info("# Using Regex patten:\n %s" % feature_regex)

    # open output string buffer (will write to compressed file later)
    reads_trimmed = StringIO.StringIO()
    reads_untrimmed = StringIO.StringIO()
    mated_reads_buffer = StringIO.StringIO()

    # Keep track of matched read IDs
    # read_ids = []

    # Find all reads containing the sequence of interest
    if input_file.endswith('.gz'):
        fastq = gzip.open(input_file, 'rb')
        fastq_mated = gzip.open(input_file_mated, 'rb')
    elif input_file.endswith('.xz'):
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
            import pdb; pdb.set_trace();

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

        # skip reads that are less than the required amount after trimming
        if len(trimmed_read[SEQUENCE_IDX]) < minimum_trimmed_length:
            continue

        # length of portion trimmed off
        trimmed_part_length = (len(read[SEQUENCE_IDX]) - 
                               len(trimmed_read[SEQUENCE_IDX]))

        # for internal matches, skip reads where match is not close enough to
        # the edge of the read
        if (trimmed_part_length - match_length) > max_dist_from_edge:
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

    log_handle.info("# Finished processing %s" % os.path.basename(input_file))

def compute_coordinates(target_genome, target_gff, feature_name, build_dir, sample_id,
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
        matched_seq         sequence matched in previous step
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
        'gene_strand', 'read_strand', 'trimmed_start', 'trimmed_stop',
        'untrimmed_seq', 'trimmed_part', 'trimmed_genome_seq', 'matched_seq',
        'matched_genome_seq', 'acceptor_site'])

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

    # Keep track of matched read IDs
    read_ids = []

    # Get coordinate and strand for each read in bam file
    for read in sam:
        # Get read where feature sequence was found
        if ((read.is_read1 and read_num == 'R2') or
            (read.is_read2 and read_num == 'R1')):
            continue

        # TESTING 2015/04/13
        if read.qname in read_ids:
            import pdb; pdb.set_trace()

        read_ids.append(read.qname)

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
                read.rlen + 1:read.pos + read.rlen + trimmed_part_length + 1].seq
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

        # Check to see that match differs from genome sequence (quick)
        if (matched_seq == matched_genome_seq):
            continue

        # Note: for real acceptor sites, actual_strand = gene_strand

        # For Poly(A) tail, extend detected read to include any A's present in
        # the genome that were trimmed off
        if feature_name == 'polya':
            match = re.search('^A*', matched_genome_seq)
            overlap_length = match.end() - match.start()

            if overlap_length > 0:
                matched_seq = matched_seq[overlap_length:]
                matched_genome_seq = matched_genome_seq[overlap_length:]

                if actual_strand == '+':
                    # give back some A's and update the relevant sequences
                    #matched_genome_seq = matched_genome_seq[overlap_length:]
                    acceptor_site = acceptor_site + overlap_length
                elif actual_strand == '-':
                    # give back some A's and update the relevant sequences
                    #matched_genome_seq = matched_genome_seq[-overlap_length:]
                    acceptor_site = acceptor_site - overlap_length

        # Poly(T)
        elif feature_name == 'polyt':
            match = re.search('T*$', matched_genome_seq)
            overlap_length = match.end() - match.start()

            if overlap_length > 0:
                matched_seq = matched_seq[overlap_length:]

                # - strand
                if actual_strand == '-':
                    # give back some A's and update the relevant sequences
                    matched_genome_seq = matched_genome_seq[:-overlap_length]
                    acceptor_site = acceptor_site - overlap_length
                # + strand (for poly(t) reads this corresponds to gene on
                # negative strand)
                elif actual_strand == '+':
                    # give back some A's and update the relevant sequences
                    matched_genome_seq = matched_genome_seq[:-overlap_length]
                    acceptor_site = acceptor_site + overlap_length

        # 
        # More filtering
        #
        # Check to see if feature is still sufficiently long;
        if len(matched_seq) < min_feature_length:
            continue

        # Check to see that match differs from genome sequence (quick)
        if (matched_seq == matched_genome_seq):
            continue

        try:
            # Make sure that for the region matched, there are at least n
            # differences between the match and the corresponding genomic
            # sequence (slower)
            seq_dist = distance.hamming(matched_seq, matched_genome_seq)
            if seq_dist < minimum_differences:
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
        if is_inside_cds(chromosomes[chromosome], acceptor_site, window_size):
            inside_cds.writerow([read.qname, chromosome, strand, acceptor_site])
            num_inside_cds = num_inside_cds + 1
            continue

        # Find nearest gene
        gene = find_closest_gene(chromosomes[chromosome], feature_name,
                                 acceptor_site, acceptor_site_side, window_size)

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
            read.qname, gene['id'], chromosome, actual_strand, strand,
            read.pos, read.pos + read.rlen, untrimmed_read.seq, trimmed_part,
            trimmed_genome_seq, matched_seq, matched_genome_seq, acceptor_site
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
            results[chromosome][gene['id']][acceptor_site]['count'] += len(matched_seq)

    # record number of good and bad reads
    log_handle.info(
        "# Found %d reads with predicted acceptor site at expected location"
        % num_good)
    log_handle.info(
        "# Found %d reads with predicted acceptor site inside a known CDS"
        % num_inside_cds)
    log_handle.info(
        "# Found %d reads with predicted acceptor site not proximal to any CDS"
        % num_no_nearby_genes)

    # Output coordinates as a GFF file
    output_filepath = '%s/%s_coordinates_%s.gff' % (
        output_dir, feature_name, read_num
    )

    output_coordinates(results, feature_name, output_filepath, target_gff)

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
    bool
        True if the feature is located in a known CDS.

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
            return True

    return False

def find_closest_gene(chromosome, feature_name, location, acceptor_site_side,
                      window_size):
    """
    Finds the closest gene to a specified location that is in the expected
    orientation.
    """
    # chromosome boundary
    ch_end =  int(chromosome.features[0].location.end)

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

        # For SL/RSL, look at gene start locations and for Poly(A)/(T) look
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
        if ((closest_gene.strand == -1 and feature_name in ['sl', 'rsl']) or 
            (closest_gene.strand ==  1 and feature_name in ['polya', 'polyt'])):
            return None
    elif acceptor_site_side == 'right':
        if ((closest_gene.strand == -1 and feature_name in ['polya', 'polyt']) or 
            (closest_gene.strand ==  1 and feature_name in ['sl', 'rsl'])):
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
