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

outfile=${1/.*}.bedgraph

# grab and comments
grep '^#' $1 > $outfile

# add track
echo 'track type=bedGraph' >> $outfile

# convert the remaining fields and append to output
grep -v '^#' $1 | awk '{print $1 "\t" ($4 - 1) "\t" $5 "\t" $6}' >> $outfile
