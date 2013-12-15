SL analysis
2013/12/12

Overview
--------

## Input

- RNA-Seq reads (fastq)
- Genome sequence (fasta)
- Genome annotation including CDS coordinates (gff)
- Spliced leader sequence

## Goal

Generate a table containing the location of primary and alternative splice
acceptor sites for each gene (and condition/sample); usage frequency for acceptor sites.

## Basic process

1. Given SL sequence, find all paired reads containing >= n bases of the SL
   sequence (current default n=10).
2. Remove reads where matching sequence is internal (SL must be at either end
   of a read)
3. Remove SL portion of sequence and map back to original genome -- site of
   mapping indicates location of splice acceptor site.

Background
----------

## RNA-Seq paired end reads

    [100bp]--(~200-400bp)--[100bp]

     between 98-100

## Gene structure (after trans-splicing)

    [SL](---5'UTR----)[CDS](3'UTR)[PolyA tail]

## SL Sequences

### T. cruzi (Yuan; source unknown)

Possible origin:
http://www.ncbi.nlm.nih.gov/nucleotide/20977244?report=genbank&log$=nuclalign&blast_rank=34&RID=ASXAKGE901R

> AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG

### L. major (Rastrojo et al. 2013)

> AACTAACGCT ATATAAGTAT CAGTTTCTGT ACTTTATTG

