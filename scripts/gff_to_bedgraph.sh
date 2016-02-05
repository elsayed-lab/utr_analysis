#!/usr/bin/env bash
#----------------------------------------------------------
# GFF to Bedgraph converter
# Keith Hughitt (khughitt@umd.edu)
# 2014/03/24
#
# Based on suggestion at:
# https://groups.google.com/forum/#!msg/igv-help/kuZk_XhbTr0/FLuuDjSd0TMJ
#
# Usage:
#  gff_to_bedgraph.sh input.gff
#----------------------------------------------------------

sorted_gff=${1/.*}_sorted.gff
outfile=${1/.*}_sorted.bedgraph

# remove ncRNAs in L. major
ncrnas=${1/.*}_no_ncrnas.gff
sorted_gff=${1/.*}_no_ncrnas_sorted.gff
grep -v "LmjF\.[0-9]*\.5\?[a-zA-Z][\.]*" ${1} > $ncrnas
igvtools sort $ncrnas $sorted_gff

#igvtools sort ${1} $sorted_gff
rm $ncrnas

# let's first sort and index the GFF files for comparison
igvtools index $sorted_gff

# remove unsorted version
rm ${1}

# grab and comments
grep '^#' $sorted_gff > $outfile

# add track
echo 'track type=bedGraph' >> $outfile

# convert the remaining fields and append to output
grep -v '^#' $sorted_gff | grep -v 'chromosome' |\
    awk '{print $1 "\t" ($4 - 1) "\t" $5 "\t" $6}' >> $outfile

rm igv.log
