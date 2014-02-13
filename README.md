UTR analysis pipeline

Overview
--------

## Input

- RNA-Seq reads (fastq)
- Genome sequence (fasta)
- Genome annotation including CDS coordinates (gff)
- Spliced leader sequence (for 5'UTR/SL analysis)

## Goals

### Spliced leader analysis

Generate a table containing the location of primary and alternative splice
acceptor sites for each gene (and condition/sample) along with the usage 
frequency for acceptor sites.

### Poly-adenylation analysis

## Basic process

1. Given SL sequence, find all paired reads containing >= n bases of the SL
   sequence (current default n=10).
2. Remove reads where matching sequence is internal (SL must be at either end
   of a read)
   - UPDATE: This may not be the case! It appears that many valid SL read
   contain unrelated sequence at the left of the read (adapter?) that make this
   filtering step overly stringent.
   - TODO: Check to see if Trimming adapter sequence improves this picture.
3. Remove SL portion of sequence and map back to original genome -- site of
   mapping indicates location of splice acceptor site.

Background
----------

## Spliced Leader

### Organization in the genome

- Tryp. genomes contain tandem arrays of several hundred SL genes from which
the SL is transcribed.
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

