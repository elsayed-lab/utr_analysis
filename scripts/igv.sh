#!/usr/bin/env bash
ANNOTATIONS=${REF}/tcruzi_clbrener/annotation/tc_esmer
READS=${HOME}/Dropbox/research/2013/20-utr-analysis/build/tophat/mismatches-0_minlength-10

${HOME}/software/IGV_2.3.26/igv.sh \
    -g ${REF}/tcruzi_clbrener/genome/tc_esmer/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like_Genome.fasta \
    ${ANNOTATIONS}/TriTrypDB-6.0_TcruziCLBrenerEsmeraldo-like.bed,${READS}/R1_match/accepted_hits_sorted.bam,${READS}/R2_match/accepted_hits_sorted.bam
