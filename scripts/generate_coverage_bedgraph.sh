#!/usr/bin/env bash
# Generates a single bedgraph representing the coverage for all of the RNA-Seq
# samples.
#
# Example usage:
# ./generate_coverage_bedgraph.sh genome.gff "/path/to/*/accepted_hits.bam"
#
# Note: wildcard string must be surround by quotations marks to be handled
# properly.

# Input arguments
GENOME_GFF=$1
INPUT_FILES=$2

# First, generate a genome file with chromosome lengths, as required by
# bedtools
grep "sequence-region" $GENOME_GFF | \
    awk -v OFS='\t' '{print $2,$4}' > /tmp/bedtools.genome

# Generate coverage map for each sample
for FILEPATH in $INPUT_FILES; do
    echo "Processing $FILEPATH"
    OUTFILE=${FILEPATH%.bam}_coverage.bedgraph

    if [ ! -e "$OUTFILE" ]; then
        bedtools genomecov -ibam -d -bga -trackline \
            -i $FILEPATH -g /tmp/bedtools.genome > $OUTFILE
    fi
done

# Combine bedgraphs into a single summary bedgraph file
echo "Combining multiple bedgraphs into a single summary bedgraph"
unionBedGraphs -empty -g /tmp/bedtools.genome \
    -i $(echo $INPUT_FILES | sed 's/\.bam/_coverage\.bedgraph/g') > coverage.bedgraph

# Generate binary TDF version as well
echo "Converting to TDF"
igvtools toTDF coverage.bedgraph coverage.tdf /tmp/bedtools.genome

