UTR Analysis Pipeline
=====================

Status: IN DEVELOPMENT

## Overview

Pipeline for analyzing a set of RNA-Seq reads to determine 5' and
3'UTR structure including the locations of all spliced leader acceptor sites
and poly-A tail acceptor sites.

The pipeline makes use of Tophat for mapping reads, and Ruffus for pipeline
flow management.

For more information on the development of this pipeline and motivation for
some of the design decisions, see my [research in progress presenation from
2013/12/20](http://www.umiacs.umd.edu/~keith/research/presentations/2013-12-20_RIP_UTR_Boundaries_Analysis/README.html).

## Input

- RNA-Seq reads (fastq)
- Genome sequence (fasta)
- Genome annotation including CDS coordinates (gff)
- Spliced leader sequence (for 5'UTR/SL analysis)
- [Optional] Genome sequen for a second species that is not of interest but may
    be used to filter out reads, e.g. in a mixed transcriptome study (fasta)
- [Optional] A GFF containing addition ORFs outside of the primary genome
    annotations to be included when assigning acceptor sites.

## Goals

### Spliced leader analysis

Generate a table containing the location of primary and alternative splice
acceptor sites for each gene (and condition/sample) along with the usage 
frequency for acceptor sites.

### Poly-adenylation analysis

## Basic process

1. Given SL sequence, find all paired reads containing >= n bases of the SL
   sequence (current default n=10).
2. Remove reads where matching sequence is internal (Optional; although in
   theory all SL-containing reads should have the SL fragment at the upstream
   end of the read, in practice there are many chimeric reads with unrelated
   fragments further upstream resulting in valid internal SL sequences.)
3. Remove SL portion of sequence and map back to original genome -- site of
   mapping indicates location of splice acceptor site.

Background
----------

## Spliced Leader

### Organization in the genome

- Trypanosome genomes contain tandem arrays of several hundred SL genes from
    which the SL is transcribed.
- "SL repeat" = SL RNA and a non-transcribed spacer.
- The SL exon sequence is highly conserved, while the spacer regions are 
  highly variable and sometimes used to differentiate closely related tryp.
  species.

### SL RNA

- The SL mini-exon is part of a larger SL RNA transcript which is thought to
fold into a 3 stem-loop structure, possibly with varying configurations
(Gibson).


### SL exon

SL = SL exon = SL mini-exon

#### T. cruzi (McCarthy-Burke et al; 1989)

> AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG

#### L. major (Rastrojo et al. 2013)

> AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG

References
----------
- Agabian, N. (1990). Trans splicing of nuclear pre-mRNAs. Cell, 61(7),
1157–1160. doi:10.1016/0092-8674(90)90674-4

- Gibson, W., Bingle, L., Blendeman, W., Brown, J., Wood, J., & Stevens, J.
(2000). Structure and sequence variation of the trypanosome spliced leader
transcript. Molecular and biochemical parasitology, 107(2), 269–77. Retrieved
from http://www.ncbi.nlm.nih.gov/pubmed/10779603

- McCarthy-Burke, C., Taylor, Z. a, & Buck, G. a. (1989). Characterization of
the spliced leader genes and transcripts in Trypanosoma cruzi. Gene, 82(1),
177–89. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/2684773

