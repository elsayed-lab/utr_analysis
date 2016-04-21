"""
File parsing and I/O functions:

    - parse_input
    - create_build_dirs
    - readfq
    - output_coordinates
    - combine_gff_results
    - load_annotations

"""
import argparse
import csv
import os
import textwrap
from BCBio import GFF

def parse_input():
    """
    Parses script input and returns values.
    """
    # Usage example
    usage_examples=textwrap.dedent("""\
    Usage Example:
    --------------
    ./utr_analysis.py                                               \\
        -i "$RAW/tcruzir21/*/processed/*.filtered.fastq.gz"         \\
        -s AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG                  \\
        -f1 TriTrypDB-27_TcruziCLBrenerEsmeraldo-like_Genome.fasta  \\
        -f2 mm10.fasta                                              \\
        -g TrypDB-27_TcruziCLBrenerEsmeraldo-like.gff               \\
        --build-directory build/tcruzi                              \\
        --min-sl-length 12                                          \\
        --exclude-internal-polya-matches
    """)

    # Create ArgumentParser instance
    parser = argparse.ArgumentParser(
        description='Spliced Leader and poly-adenylation site analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=usage_examples
    )

    # Add arguments
    parser.add_argument('-i', '--input-reads', required=True,
                        help=('RNA-Seq FASTQ or gzipped/xz-compressed FASTQ '
                              'glob string or a txt file containing filepaths '
                              'to the samples to be used'))
    parser.add_argument('-d', '--build-directory', required=True,
                        help='Directory to save output to')
    parser.add_argument('-f1', '--target-genome', dest='target_genome',
                        required=True, help=('Genome sequence FASTA filepath '
                        'for target species'))
    parser.add_argument('-f2', '--host-genome', dest='host_genome',
                        help=('Genome sequence FASTA filepath for species to '
                              'be filtered out prior to mapping. (optional)'))
    parser.add_argument('-g1', '--target-annotations', dest='target_gff',
                        required=True, help='Genome annotation GFF')
    parser.add_argument('-g2', '--host-annotations', dest='host_gff',
                        help='Genome annotation GFF')
    parser.add_argument('-s', '--sl-sequence', dest='spliced_leader',
                        required=True, help='Spliced leader DNA sequence',
                        default=None)
    parser.add_argument('--exclude-internal-sl-matches',
                        help=('Only allow matches with the SL at the upstream '
                              'end of a read.'), action='store_true')
    parser.add_argument('--exclude-internal-polya-matches',
                        help=('Only allow matches with the Poly(A) tail at the '
                              'downstream end of a read.'), action='store_true')
    parser.add_argument('--minimum-trimmed-length',
                        help=('The minimum read length allowed after SL/Poly(A)'
                              'trimming has been performed. (default=20)'),
                        default=20)
    parser.add_argument('--max-dist-from-edge',
                        help=('For unanchored searches, what is the maximum '
                        'distance from the edge of the read for a feature '
                        'match to be considered (default=unlimited).'),
                        default=[])
    parser.add_argument('-m', '--min-sl-length', default=6, type=int,
                        help='Minimum length of SL match (default=6)')
    parser.add_argument('-p', '--min-polya-length', default=6, type=int,
                        help='Minimum length of Poly-A match (default=6)')
    parser.add_argument('-w', '--window-size', default=50000, type=int,
                        help=('Number of bases up or downstream of feature to '
                              'scan for related genes (default=50000)'))
    parser.add_argument('-x', '--minimum-differences', default=2, type=int,
                        help=('Minimum number of differences from genomic '
                              ' sequence for a hit to be considered real. '
                              '(default=2)'))
    parser.add_argument('--num-threads', default=4, type=int,
                        help='Number of threads to use (default=4).')
    parser.add_argument('--num-threads-tophat', default=1, type=int,
                        help=('Number of threads to use for each Tophat run. '
                              '(default=1)'))

    # Parse arguments
    args = parser.parse_args()

    # Replace any environmental variables and return args
    args.target_genome = os.path.expandvars(args.target_genome)
    args.input_reads = os.path.expandvars(args.input_reads)
    args.target_gff = os.path.expandvars(args.target_gff)

    if args.host_genome:
        args.host_genome = os.path.expandvars(args.host_genome)
    if args.host_gff:
        args.host_gff = os.path.expandvars(args.host_gff)

    # @TODO Validate input
    return args

def create_build_dirs(args, sample_ids):
    """Creates necessary build directories for the pipeline and returns
       a dictionary of the directory paths"""
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
    sl_dir_suffix = 'minlength%d%s_mindiff-%d' % (
            args.min_sl_length,
            '-anchored' if args.exclude_internal_sl_matches else '',
            args.minimum_differences)

    # SL sub-directory
    sl_build_dir = os.path.join(args.build_directory,
                                'spliced_leader',
                                sl_dir_suffix)

    # Reverse SL sub-directory
    rsl_build_dir = os.path.join(args.build_directory,
                                    'reverse_spliced_leader',
                                    sl_dir_suffix)

    # Poly(A) tail sub-directory
    polya_dir_suffix = 'minlength%d%s_mindiff-%d' % (
            args.min_polya_length,
            '-anchored' if args.exclude_internal_polya_matches else '',
            args.minimum_differences)

    polya_build_dir = os.path.join(args.build_directory,
                                    'poly-a',
                                    polya_dir_suffix)

    # Poly(A) tail reverse complement sub-directory
    polyt_build_dir = os.path.join(args.build_directory,
                                   'poly-t',
                                    polya_dir_suffix)

    # create subdirs based on matching parameters
    shared_subdirs = ['fastq', 'tophat/mapped_to_target_untrimmed']

    if args.host_genome:
        shared_subdirs = shared_subdirs + ['tophat/mapped_to_host']

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
                            'results', 'log', 'tophat']:
                outdir = os.path.join(base_dir, sample_id, sub_dir)
                if not os.path.exists(outdir):
                    os.makedirs(outdir, mode=0o755)

    return {
        'shared': shared_build_dir,
        'sl': sl_build_dir,
        'rsl': rsl_build_dir,
        'polya': polya_build_dir,
        'polyt': polyt_build_dir,
        'combined': combined_output_dir
    }

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

def output_coordinates(results, feature_name, filepath, target_gff, track_color='0,0,255'):
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
    annotations_fp = open(target_gff)

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
                attributes = "ID=%s.%s.%d;Name=%s;seq_dist=%d;match_len=%d;description=%s" % (
                    gene_id, feature_name, i, gene_id,
                    results[chrnum][gene_id][acceptor_site]['seq_dist'],
                    results[chrnum][gene_id][acceptor_site]['match_len'],
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
            # Example:
            # 'ID=TcCLB.507939.50.sl.2;Name=TcCLB.507939.50;seq_dist=10;
            #  match_len=15;description=hypothetical+protein,+conserved'
            parts = row['description'].split(';')

            NAME_IDX = 1
            SEQDIST_IDX = 2
            MATCHLEN_IDX = 3
            DESC_IDX = 4

            # Parse individual attributes

            gene = parts[NAME_IDX].split('=').pop()
            seq_dist = int(parts[SEQDIST_IDX].split('=').pop())
            match_len = int(parts[MATCHLEN_IDX].split('=').pop())
            description = parts[DESC_IDX].split('=').pop()

            # Add to output dictionary
            if not chromosome in results:
                results[chromosome] = {}
            if not gene in results[chromosome]:
                results[chromosome][gene] = {}

            # Increment site count and save distance from gene
            if not acceptor_site in results[chromosome][gene]:
                results[chromosome][gene][acceptor_site] = {
                    "count": count,
                    "seq_dist": seq_dist,
                    "match_len": match_len,
                    "strand": strand,
                    "description": description
                }
            else:
                results[chromosome][gene][acceptor_site]['count'] += count
                results[chromosome][gene][acceptor_site]['seq_dist'] += seq_dist
                results[chromosome][gene][acceptor_site]['match_len'] += match_len

    # Return combined results
    return results

def load_annotations(target_gff):
    """Loads genome annotations from specified GFF(s)."""
    # Get chromosomes/contigs from GFF file
    chromosomes = {}

    # Load existing gene annotations
    annotations_fp = open(target_gff)

    for entry in GFF.parse(annotations_fp):
        if len(entry.features) > 0 and entry.features[0].type in ['chromosome', 'contig']:
            chromosomes[entry.id] = entry

    # clean up
    annotations_fp.close()

    return chromosomes

