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
        fp.writelines(utr5_lengths_proc + utr5_lengths_meta)

    # 3'UTR lengths
    printb("Computing 3'UTR lengths...")

    utr3_lengths_proc = compute_3utr_lengths(polya_proc, genes,
                                         'output/lmajor_3UTR_lengths_proc.txt')
    utr3_lengths_meta = compute_3utr_lengths(polya_meta, genes,
                                         'output/lmajor_3UTR_lengths_meta.txt')
    # Write combined 3'UTR lengths
    with open('output/lmajor_3UTR_lengths_all.txt', 'w') as fp:
        fp.writelines(utr3_lengths_proc + utr3_lengths_meta)

    # Primary/minor SL acceptor site distances
    printb("Computing distances between primary and minor SL acceptor sites...")

    sl_site_distances_proc = compute_alt_acceptor_site_distances(
        sl_proc, genes, 'output/lmajor_dist_sl_sites_proc.txt')
    sl_site_distances_meta = compute_alt_acceptor_site_distances(
        sl_meta, genes, 'output/lmajor_dist_sl_sites_meta.txt')
    with open('output/lmajor_dist_sl_sites_all.txt', 'w') as fp:
        fp.writelines(sl_site_distances_proc + sl_site_distances_meta)

    # Primary/minor Poly(A) acceptor site distances
    printb("Computing distances between primary and minor Poly(A) acceptor sites...")

    polya_site_distances_proc = compute_alt_acceptor_site_distances(
        polya_proc, genes, 'output/lmajor_dist_polya_sites_proc.txt')
    polya_site_distances_meta = compute_alt_acceptor_site_distances(
        polya_meta, genes, 'output/lmajor_dist_polya_sites_meta.txt')
    with open('output/lmajor_dist_polya_sites_all.txt', 'w') as fp:
        fp.writelines(polya_site_distances_proc + polya_site_distances_meta)

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

            # otherwise get the UTR length
            if gene.strand == 1:
                length = gene.location.start - entry.location.start
            else:
                length = entry.location.start - gene.location.end

            utr5_lengths.append("%s\n" % length)

    # save output
    with open(outfile, 'w') as fp:
        fp.writelines(utr5_lengths)

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

            # otherwise get the UTR length
            if gene.strand == 1:
                length = entry.location.start - gene.location.end
            else:
                length = gene.location.start - entry.location.start

            utr3_lengths.append("%s\n" % length)

    # save output
    with open(outfile, 'w') as fp:
        fp.writelines(utr3_lengths)

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
        gene_ids = list(set(gene_ids))

        # iterate over genes
        for gene_id in gene_ids:
            sites = [x for x in ch.features if x.qualifiers['Name'][0] ==
                     gene_id]

            # if only primary site found, stop here
            if len(sites) == 1:
                continue

            # otherwise, determine primary site
            primary = sites[0]
            max_coverage = int(primary.qualifiers['score'][0])

            for site in sites[1:]:
                coverage = int(site.qualifiers['score'][0])
                if coverage > max_coverage:
                    primary = site
                    max_coverage = coverage

            # Measure distances from primary site
            for site in sites:
                strand = genes[gene_id].strand

                if strand == 1:
                    dist = site.location.start - primary.location.start
                else:
                    dist = primary.location.start - site.location.start

                if dist != 0:
                    distances.append("%s\n" % dist)

    with open(outfile, 'w') as fp:
        fp.writelines(distances)

    return distances

"""
MAIN
"""
if __name__ == "__main__":
    main()

