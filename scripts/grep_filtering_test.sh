#!/usr/bin/env bash
#
# Performance test: FASTQ filtering with grep
# Keith
# 2014/02/11
#
# The goal of this script is to test the performance of using grep to filter
# a large (e.g. 10s of millions of records) FASTQ file.
# Currently this is done in Python and is one of the bottle-necks in the
# processing pipeline.

# Settings
N=135000
INPUT=$RAW/tcruzir21/HPGL0258/processed/HPGL0258_R1_combined.filtered.fastq
INTERMEDIATE="random_ids.$N"
OUTPUT=grep_test.fastq

# Generate a random subset of the reads
if [ ! -e "$INTERMEDIATE" ]; then
    echo "Generating random subset of reads..."
    grep --color='never' HWI $INPUT | \
        awk '{print $1}' | cut -c '2-' | sort -R | head -n$N > $INTERMEDIATE
fi

COUNTER=1

if [ -e $OUTPUT ]; then
    rm $OUTPUT;
fi

echo "Grepping for those ids in original file"
t1=$(date +%s)
while read line; do
    printf "%d/%d\n" $COUNTER $N
    grep --color='never' -A 3 $line $INPUT >> $OUTPUT
    COUNTER=$[$COUNTER +1]
done < $INTERMEDIATE
t2=$(date +%s)

echo "Done!"
echo "Finished in " $(expr $t2 - $t1) "Seconds..."

