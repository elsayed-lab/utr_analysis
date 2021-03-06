#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
Intergenic ORF Detection
Keith Hughitt (khughitt@umd.edu)
2014/03/25

A simple script to look for ORFs of some specified minimum length in between 
known CDSs.

Usage Example: 

    ./orf_detection.py genome.fasta annotations.gff 50

Finds all ORFs of at least 50 peptides that fall outside of known CDSs in the 
specified genome and annotations files.
"""
import os
import csv
import sys
from Bio import Seq, SeqIO
from BCBio import GFF

def main():
    """Main"""
    # Parse input
    (fasta, gff, min_length) = parse_input()

    # Load genome
    chromosomes = {x.id:x for x in SeqIO.parse(fasta, 'fasta')}

    # Load annotations
    annotations = {}

    for entry in GFF.parse(open(gff)):
        if len(entry.features) > 0 and entry.features[0].type == 'chromosome':
            annotations[entry.id] = entry

    # Determine locations of inter-CDS regions for each chromosome
    inter_cds_regions = {}

    for chrnum, chromosome in annotations.items():
        # get chromosome dimensions (the first feature represents the
        # chromosome itself)
        ch_end = int(chromosome.features[0].location.end)

        # filter out everything except for genes
        genes = [x for x in chromosome.features if x.type == 'gene']

        # add range before first gene
        start = 0

        # keep track of strand of polycistronic transcriptional unit
        strand = None

        inter_cds_regions[chrnum] = {
            -1: [],
            +1: []
        }

        # iterate through genes and store the ranges between them;
        # for TriTrypDB files, the gene boundaries are the same as the CDS
        # boundaries (@TODO: check if this is also true for T. brucei)
        for gene in genes:
            end = int(gene.location.start)

            # Add CDS to relevant list based on strand
            if strand is None:
                # Left-most gene
                inter_cds_regions[chrnum][gene.location.strand].append((start, end))
            elif strand != gene.location.strand:
                # Add ORFs in both directions at transcription switch sites (TSSs)
                inter_cds_regions[chrnum][+1].append((start, end))
                inter_cds_regions[chrnum][-1].append((start, end))
            else:
                # Within PTU
                inter_cds_regions[chrnum][strand].append((start, end))

            # update start counter and strand
            start = int(gene.location.end)
            strand = gene.location.strand

        # add region after last gene
        inter_cds_regions[chrnum][strand].append((start, ch_end))

    # Iterate over inter-CDS regions and find ORFs of at least the specified 
    # length in any of the six possible reading frames and output as GFF entries
    fp = open('output/orfs.gff', 'w')

    # Write csv header
    fp.write("##gff-version\t3\n")
    fp.write("##feature-ontology\tsofa.obo\n")
    fp.write("##attribute-ontology\tgff3_attributes.obo\n")

    # Write header to output
    writer = csv.writer(fp, delimiter='\t')

    # Iterate through interCDS regions
    for i, chrnum in enumerate(inter_cds_regions, 1):
        print("Processing %s (%d/%d)" % (chrnum, i, len(chromosomes)))
        for strand in [-1, 1]:
            for j, region in enumerate(inter_cds_regions[chrnum][strand]):
                # get sequence record for the range
                record = chromosomes[chrnum][region[0]:region[1]]

                # Find ORFs in each of the six possible reading frames 
                # that are at least the specified length
                orfs = find_orfs(record.seq, min_length, strand)

                # Write GFF entries for each match
                for k, orf in enumerate(orfs, 1):
                    # Convert coordinates to chromosomal position
                    #if orf[2] == "+":
                    start = orf[0] + region[0] + 1
                    stop = orf[1] + region[0]
                    str_strand = orf[2]

                    # Write entry
                    orf_id = "%s.ORF.%d.%d" % (chrnum, j, k)
                    gff_attrs = "ID=%s;Name=%s" % (orf_id, orf_id)
                    writer.writerow([chrnum, "ElSayedLab", 'ORF',
                                start, stop, '.', str_strand, '.', gff_attrs])

    # clean up
    fp.close()

def parse_input():
    """Parses and validates input"""
    # check number of args
    if not len(sys.argv) == 4:
        print("Incorrect number of arguments specified!")
        print("Example usage:")
        print("orf_detection.py genome.fasta annotations.gff 50")
        sys.exit(0)

    (fasta, gff, min_length) = sys.argv[1:]

    # case minimum length to an integer
    min_length = int(min_length)

    # check to make sure valid filepaths specified
    if not os.path.exists(fasta):
        print("Incorrect genome filepath specified")
        sys.exit(0)
    if not os.path.exists(gff):
        print("Incorrect annotations filepath specified")
        sys.exit(0)

    # return input arguments
    return (fasta, gff, min_length)

    def find_orfs(seq, min_protein_length, strand=1, trans_table=1,
                  ignore_ambiguous_orfs=True):
        """
        Finds ORFs of a specified minimum protein length in a SeqRecord.

        Based on: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec360
        """
        answer = []
        seq_len = len(seq)

        # Get sequence associated with the specified location and strand
        if strand == 1:
            dna_seq = seq
        else:
            dna_seq = seq.reverse_complement()

        for frame in range(3):
            trans = str(dna_seq[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0

            # Iterate through ORFS in reading frame
            while aa_start < trans_len:
                # Set end counter to position of next stop codon
                aa_start = trans.find("M", aa_start)
                aa_end = trans.find("*", aa_start)

                # If no start or stop codons found, stop here
                if aa_start == -1 or aa_end == -1:
                    break

                if (aa_end < aa_start):
                    raise Exception('wtf')

                # Compute coordinates of ORF
                if strand == 1:
                    start = frame + aa_start * 3
                    end = min(seq_len, frame + aa_end * 3 + 3)
                else:
                    start = seq_len - frame - aa_end * 3 - 3
                    end = seq_len - frame - aa_start * 3

                # Add to output
                str_strand = "+" if strand == 1 else '-'

                # Check to make sure ORF doesn't contain a bunch of N's
                if ignore_ambiguous_orfs:
                    num_unknown = trans[aa_start:aa_end].count('X')
                    if (num_unknown / (aa_end - aa_start)) > 0.25:
                        aa_start = aa_end + 1
                        continue

                # increment start counter
                aa_start = aa_end + 1

                # Add ORF coordinates and continue 
                answer.append((start, end, str_strand))

        # Sort results
        answer.sort()

        return answer

if __name__ == "__main__":
    main()
