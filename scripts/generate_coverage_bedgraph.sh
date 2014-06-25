#!/usr/bin/env bash
# Generates a single bedgraph representing the coverage for all of the RNA-Seq
# samples.
#
# Example usage:
# ./generate_coverage_bedgraph.sh genome.gff "/path/to/*/accepted_hits.bam"
# ./generate_coverage_bedgraph.sh hg19 "/path/to/*/accepted_hits.bam"
#
# Note: wildcard string must be surround by quotations marks to be handled
# properly.
#
# Note 2: Code may be adjusted to use actual human.genome file instead of
# just 'hg19' for the first couple steps when generating coverage file for
# human.

# Input arguments
GENOME=$1
INPUT_FILES=$2

# First, generate a genome file with chromosome lengths, as required by
# bedtools
if [ "${GENOME##*.}" = "gff" ]; then
    GENOME_SIZES=/tmp/bedtools.genome
    grep "sequence-region" $GENOME | \
        awk -v OFS='\t' '{print $2,$4}' > $GENOME_SIZES
else
    GENOME_SIZES=$GENOME
fi

# Generate coverage map for each sample
for FILEPATH in $INPUT_FILES; do
    echo "Processing $FILEPATH"
    OUTFILE=${FILEPATH%.bam}_coverage.bedgraph

    if [ ! -e "$OUTFILE" ]; then
        bedtools genomecov -ibam -d -bga -trackline \
            -i $FILEPATH -g $GENOME_SIZES > $OUTFILE
    fi
done

# Combine bedgraphs into a single summary bedgraph file
echo "Combining multiple bedgraphs into a single summary bedgraph"
unionBedGraphs -empty -g $GENOME_SIZES \
    -i $(echo $INPUT_FILES | sed 's/\.bam/_coverage\.bedgraph/g') > coverage.bedgraph

# Delete individual bedgraph files

# Generate binary TDF version as well
echo "Converting to TDF"
igvtools toTDF coverage.bedgraph coverage.tdf $GENOME_SIZES

