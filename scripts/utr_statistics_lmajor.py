#!/usr/bin/env python
"""
UTR statistics
2014/06/27
Keith Hughitt (khughitt@umd.edu)

Parses results from UTR analysis pipeline and generates some derrived data
files and statistics on UTR length, acceptor site usage, etc. to be used in
downstream analysis scripts by YL.

Usage
-----
./utr_statistics.py genome.gff spliced_leader.gff polya.gff
"""
import os
import csv
import sys
from BCBio import GFF

def main():
    """Main"""
    # Procyclic
    procyclic_dir = os.path.expanduser(
        '~/sl_analysis_build/lmajor-procyclic/results/minlength-4/unanchored')

    sl_gff_proc = os.path.join(procyclic_dir, 'spliced_leader.gff')
    polya_gff_proc = os.path.join(procyclic_dir, 'polya.gff')

    # Metacyclic
    metacyclic_dir = os.path.expanduser(
        '~/sl_analysis_build/lmajor-metacyclic/results/minlength-4/unanchored')

    sl_gff_meta = os.path.join(metacyclic_dir, 'spliced_leader.gff')
    polya_gff_meta = os.path.join(metacyclic_dir, 'polya.gff')

    # Load gene annotations
    printb("Loading Genome annotations...")

    gff = "/cbcb/lab/nelsayed/ref_data/lmajor_friedlin/annotation/TriTrypDB-8.0_LmajorFriedlin.gff"

    genes = {}
    for ch in GFF.parse(open(gff)):
        for entry in ch.features:
            if entry.type == 'gene':
                genes[entry.id] = entry

    # Load SL and Poly(A) sites
    printb("Loading SL acceptor sites...")
    sl_proc = list(GFF.parse(open(sl_gff_proc)))
    sl_meta = list(GFF.parse(open(sl_gff_meta)))

    printb("Loading Poly(A) acceptor sites...")
    polya_proc = list(GFF.parse(open(polya_gff_proc)))
    polya_meta = list(GFF.parse(open(polya_gff_meta)))

    # 5'UTR lengths
    printb("Computing 5'UTR lengths...")

    utr5_lengths_proc = compute_5utr_lengths(sl_proc, genes,
                                         'output/lmajor_5UTR_lengths_proc.txt')
    utr5_lengths_meta = compute_5utr_lengths(sl_meta, genes,
                                         'output/lmajor_5UTR_lengths_meta.txt')
    # Write combined 5'UTR lengths
    with open('output/lmajor_5UTR_lengths_all.txt', 'w') as fp:
        fp.writelines(sorted(utr5_lengths_proc + utr5_lengths_meta))

    ########################################
    # TESTING
    ########################################
    sys.exit()

    # 3'UTR lengths
    printb("Computing 3'UTR lengths...")

    utr3_lengths_proc = compute_3utr_lengths(polya_proc, genes,
                                         'output/lmajor_3UTR_lengths_proc.txt')
    utr3_lengths_meta = compute_3utr_lengths(polya_meta, genes,
                                         'output/lmajor_3UTR_lengths_meta.txt')
    # Write combined 3'UTR lengths
    with open('output/lmajor_3UTR_lengths_all.txt', 'w') as fp:
        fp.writelines(sorted(utr3_lengths_proc + utr3_lengths_meta))

    # Primary/minor SL acceptor site distances
    printb("Computing distances between primary and minor SL acceptor sites...")

    sl_site_distances_proc = compute_alt_acceptor_site_distances(
        sl_proc, genes, 'output/lmajor_dist_sl_sites_proc.txt')
    sl_site_distances_meta = compute_alt_acceptor_site_distances(
        sl_meta, genes, 'output/lmajor_dist_sl_sites_meta.txt')
    with open('output/lmajor_dist_sl_sites_all.txt', 'w') as fp:
        fp.writelines(sorted(sl_site_distances_proc + sl_site_distances_meta))

    # Primary/minor Poly(A) acceptor site distances
    printb("Computing distances between primary and minor Poly(A) acceptor sites...")

    polya_site_distances_proc = compute_alt_acceptor_site_distances(
        polya_proc, genes, 'output/lmajor_dist_polya_sites_proc.txt')
    polya_site_distances_meta = compute_alt_acceptor_site_distances(
        polya_meta, genes, 'output/lmajor_dist_polya_sites_meta.txt')
    with open('output/lmajor_dist_polya_sites_all.txt', 'w') as fp:
        fp.writelines(sorted(polya_site_distances_proc +
                             polya_site_distances_meta))

    # Primary/minor SL site usage stats in procyclic and metacyclic samples
    printb("Determining primary and minor trans-splicing site usage across conditions")
    output_site_usage(genes,
                      sl_proc, 'procyclic_num_reads',
                      sl_meta, 'metacyclic_num_reads',
                      'output/lmajor_primary_trans_splicing_site_proc.txt',
                      'output/lmajor_minor_trans_splicing_site_proc.txt')

    output_site_usage(genes,
                      sl_meta, 'metacyclic_num_reads',
                      sl_proc, 'procyclic_num_reads',
                      'output/lmajor_primary_trans_splicing_site_meta.txt',
                      'output/lmajor_minor_trans_splicing_site_meta.txt')

    # Primary/minor Poly(A) site usage stats in procyclic and metacyclic samples
    printb("Determining primary and minor poly-adenylation acceptor site usage across conditions")
    output_site_usage(genes,
                      polya_proc, 'procyclic_num_reads',
                      polya_meta, 'metacyclic_num_reads',
                      'output/lmajor_primary_polya_site_proc.txt',
                      'output/lmajor_minor_polya_site_proc.txt')

    output_site_usage(genes,
                      polya_meta, 'metacyclic_num_reads',
                      polya_proc, 'procyclic_num_reads',
                      'output/lmajor_primary_polya_site_meta.txt',
                      'output/lmajor_minor_polya_site_meta.txt')

def printb(text):
    """Print text in bold"""
    print("\033[1m%s\033[0m" % text)

#
# Question 1: what is the distribution of 5' and 3'UTR lengths?
#
def compute_5utr_lengths(sl, genes, outfile):
    """Determines the UTR lengths corresponding to a set of SL sites"""
    utr5_lengths = []

    # Iterate over SL acceptor sites
    for ch in sl:
        for entry in ch.features:
            # skip chromosomes
            if entry.type == 'chromosome':
                continue

            # gene id
            gene_id = entry.qualifiers['Name'][0]

            # for now, only consider UTRs of known genes
            if gene_id.startswith('ORF'):
                continue
            gene = genes[gene_id]

            # For information on biopython feature locations, see:
            # http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html
            acceptor_site_location = entry.location.end

            # otherwise get the UTR length
            if gene.strand == 1:
                length = (gene.location.start + 1) - acceptor_site_location
            else:
                length = acceptor_site_location - gene.location.end

            if length < 1:
                import pdb; pdb.set_trace()

            utr5_lengths.append("%s\n" % length)

    # save output
    with open(outfile, 'w') as fp:
        fp.writelines(sorted(utr5_lengths))

    return utr5_lengths

def compute_3utr_lengths(polya, genes, outfile):
    """Determines the UTR lengths corresponding to a set of Poly(A) sites"""
    utr3_lengths = []

    # Iterate over Poly(A) acceptor sites
    for ch in polya:
        for entry in ch.features:
            # skip chromosomes
            if entry.type == 'chromosome':
                continue

            # gene id
            gene_id = entry.qualifiers['Name'][0]

            # for now, only consider UTRs of known genes
            if gene_id.startswith('ORF'):
                continue
            gene = genes[gene_id]

            # For information on biopython feature locations, see:
            # http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html
            acceptor_site_location = entry.location.end

            # otherwise get the UTR length
            if gene.strand == 1:
                length = acceptor_site_location - gene.location.end
            else:
                length = (gene.location.start + 1) - acceptor_site_location

            utr3_lengths.append("%s\n" % length)

    # save output
    with open(outfile, 'w') as fp:
        fp.writelines(sorted(utr3_lengths))

    return utr3_lengths

#
# Question 2: what is the distribution of distances between primary and
# alternative SL/Poly(A) acceptor sites?
#
def compute_alt_acceptor_site_distances(acceptor_sites, genes, outfile):
    """Determines the distribution of distances betwene primary and alternative
    acceptor sites"""
    # iterate over chromosomes
    distances = []

    for ch in acceptor_sites:
        # genes for which acceptor sites have been observed
        gene_ids = get_covered_geneids(ch)

        # iterate over genes
        for gene_id in gene_ids:
            # get acceptor sites for current gene
            sites = get_gene_acceptor_sites(ch, gene_id)

            # if only primary site found, stop here
            if len(sites) == 1:
                continue

            # otherwise, determine primary and minor sites
            (primary_site, minor_sites) = get_primary_and_minor_acceptor_site(sites)

            # Measure distances from primary site
            for site in minor_sites:
                strand = genes[gene_id].strand

                if strand == 1:
                    dist = site.location.start - primary_site.location.start
                else:
                    dist = primary_site.location.start - site.location.start

                distances.append("%s\n" % dist)

    with open(outfile, 'w') as fp:
        fp.writelines(sorted(distances))

    return distances

def get_gene_acceptor_sites(ch, gene_id):
    """Returns the acceptor sites specific for a gene"""
    return [x for x in ch.features if x.qualifiers['Name'][0] == gene_id]

def get_primary_and_minor_acceptor_site(sites):
    # start by assuming first site is the primary one
    primary_site = sites[0]
    primary_index = 0
    max_coverage = int(primary_site.qualifiers['score'][0])

    # iterate over remaining sites and find primary site
    for i, site in enumerate(sites[1:], 1):
        coverage = int(site.qualifiers['score'][0])
        if coverage > max_coverage:
            primary_site = site
            primary_index = i
            max_coverage = coverage

    # minor sites are all those except for the primary
    minor_sites = sites[:primary_index] + sites[(primary_index + 1):]

    return (primary_site, minor_sites)

def get_covered_geneids(ch):
    """Gets a list of gene IDs for which acceptor sites have been observed"""
    # first, get a list of all of the genes covered
    gene_ids = list(set([x.qualifiers['Name'][0] for x in ch.features
                         if x.type != 'chromosome']))
    gene_ids = []
    for entry in ch.features:
        if entry.type == 'chromosome':
            continue
        gene_id = entry.qualifiers['Name'][0]

        if gene_id.startswith('ORF'):
            continue
        gene_ids.append(gene_id)

    # unique ids
    return list(set(gene_ids))

def get_site_by_location(sites, loc):
    """Looks for an acceptor site by location and return it if found"""
    for site in sites:
        if site.location.start == loc:
            return site
    return None
#
# Question 3: How does site usage differ between procyclic and metacyclic
#
    # lmajor_primary_trans_splicing_site_proc.v2.txt
    #
    # gene              strand   chr     cds.end   acceptor.site  # proc  # meta
    # LmjF.01.0010       -       01       4703       5054       174       384
    # LmjF.01.0020       -       01       7440       7631       71       108
    # LmjF.01.0030       -       01       11068       11655       148       166

def output_site_usage(genes, conda, conda_name, condb, condb_name,
                      primary_outfile, minor_outfile):
    """Computes acceptor site usage statistics for primary and minor sites
    under varying conditions"""
    # primary site CSV writer
    writer_primary_conda = csv.writer(open(primary_outfile, 'w'),
                                      delimiter='\t')

    writer_primary_conda.writerow([
        'gene', 'strand', 'chromosome',
        'utr_start', 'utr_stop', conda_name, condb_name
    ])

    # minor site CSV writer
    writer_minor_conda = csv.writer(open(minor_outfile, 'w'), delimiter='\t')
    writer_minor_conda.writerow([
        'gene', 'strand', 'chromosome',
        'utr_start', 'utr_stop', conda_name, condb_name
    ])

    # iterate over condition_a and condition_b acceptor sites
    for i in range(len(conda)):
        ch_conda = conda[i]
        ch_condb = condb[i]

        # chromosome number
        ch_num = ch_conda.id[-2:]

        # genes for which acceptor sites have been observed (condaylic)
        gene_ids_conda = get_covered_geneids(ch_conda)

        # genes for which acceptor sites have been observed (condition_b)
        gene_ids_condb = get_covered_geneids(ch_condb)

        for gene_id in gene_ids_conda:
            # if gene only covered in one condition, stop here
            if gene_id not in gene_ids_condb:
                continue

            # get acceptor sites for current gene under each condition
            sites_conda = get_gene_acceptor_sites(ch_conda, gene_id)
            sites_condb = get_gene_acceptor_sites(ch_condb, gene_id)

            # if only primary site found, stop here
            if len(sites_conda) == 1:
                continue

            # gene strand
            gene = genes[gene_id]
            strand = gene.strand

            # otherwise, determine primary and minor sites
            (primary_site, minor_sites) = (
                get_primary_and_minor_acceptor_site(sites_conda))

            #
            # Primary site
            #

            # number of reads for condition_a primary acceptor site in condition_a
            conda_primary_site_reads_conda = primary_site.qualifiers['score'][0]

            # find same acceptor site in condition_b if it exists
            primary_site_in_condb = get_site_by_location(sites_condb,
                    primary_site.location.start)

            # condaylic primary site coverae for condition_b samples
            if primary_site_in_condb is None:
                conda_primary_site_reads_condb = 0
            else:
                conda_primary_site_reads_condb = (
                    primary_site_in_condb.qualifiers['score'][0])

            # determine UTR start and stop locations
            if strand == 1:
                start = primary_site.location.start
                stop = gene.location.start
            else:
                start = gene.location.end
                stop = primary_site.location.start

            # write row
            # LmjF.01.0010       -       01       4703       5054       174       384
            writer_primary_conda.writerow([
                gene_id, strand, ch_num, start, stop,
                conda_primary_site_reads_conda, conda_primary_site_reads_condb
            ])

            #
            # Minor sites
            #
            for minor_site in minor_sites:
                # number of reads for condition_a minor acceptor site in condition_a
                conda_minor_site_reads_conda = minor_site.qualifiers['score'][0]

                # find same acceptor site in condition_b if it exists
                minor_site_in_condb = get_site_by_location(sites_condb,
                        minor_site.location.start)

                # condaylic minor site coverae for condition_b samples
                if minor_site_in_condb is None:
                    conda_minor_site_reads_condb = 0
                else:
                    conda_minor_site_reads_condb = (
                        minor_site_in_condb.qualifiers['score'][0])

                # determine UTR start and stop locations
                if strand == 1:
                    start = minor_site.location.start
                    stop = gene.location.start
                else:
                    start = gene.location.end
                    stop = minor_site.location.start

                # write row
                writer_minor_conda.writerow([
                    gene_id, strand, ch_num, start, stop,
                    conda_minor_site_reads_conda, conda_minor_site_reads_condb
                ])

"""
MAIN
"""
if __name__ == "__main__":
    main()

