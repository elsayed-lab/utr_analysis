#!/usr/bin/env python
"""
Bedgraph subtraction
Keith Hughitt (khughitt@umd.edu)
2014/07/23

Takes two bedgraph files covering the same region (e.g. an entire genome) and
subtracts one from the other generating a "difference" bedgraph file with
positive and negative values.

This script assumes that the bedgraph files contain 5 metadata lines at the
beginning of the file which are skipped. If the input files are formatted
differently, the code can be easily adjusted below.

Usage
-----
python subtract_bedgraphs.py 1.bedgraph 2.bedgraph

To sort the output, you can then use the command:
    sort -k1,1 -k2n difference.bedgraph > difference_sorted.bedgraph
"""
import sys
import csv

# Min coverage
# Set this above one to ignore any sites with coverage below a certain value
MIN_COVERAGE = 2

# Load files
if len(sys.argv) != 3:
    print("Invalid usage: please specify two bedgraph files")
    sys.exit()

infile1, infile2 = sys.argv[1:3]

fp1 = open(infile1)
fp2 = open(infile2)

# Skip first few lines (comment and track info) for each input bedgraph
NUM_SKIP_LINES = 5

for fp in [fp1, fp2]:
    for i in range(NUM_SKIP_LINES):
        fp.readline()

# Create CSV reader instances
fields = ['chr', 'start', 'end', 'coverage']

csv1 = csv.DictReader(fp1, fieldnames=fields, delimiter='\t')
csv2 = csv.DictReader(fp2, fieldnames=fields, delimiter='\t')

# Iterate over entries in the first bedgraph and check for corresponding
# locations in the second file
difference = {}

# first bedgraph
for row in csv1:
    # skip low coverage positions
    if int(row['coverage']) < MIN_COVERAGE:
        continue

    # check for existing entry at location
    if row['chr'] not in difference:
        difference[row['chr']] = {}
    if row['start'] not in difference[row['chr']]:
        difference[row['chr']][row['start']] = {}
    if row['end'] not in difference[row['chr']][row['start']]:
        difference[row['chr']][row['start']][row['end']] = int(row['coverage'])

# second bedgraph
for row in csv2:
    # skip low coverage positions
    if int(row['coverage']) < MIN_COVERAGE:
        continue

    # check for existing entry at location
    if row['chr'] not in difference:
        difference[row['chr']] = {}
    if row['start'] not in difference[row['chr']]:
        difference[row['chr']][row['start']] = {}
    if row['end'] not in difference[row['chr']][row['start']]:
        difference[row['chr']][row['start']][row['end']] = -int(row['coverage'])
    else:
        difference[row['chr']][row['start']][row['end']] = (
            difference[row['chr']][row['start']][row['end']] -
            int(row['coverage'])
        )

# Save output
outfile = open("output/difference.bedgraph", "w")

for ch in difference:
    for start in difference[ch]:
        for end in difference[ch][start]:
            coverage = str(difference[ch][start][end])
            line = "\t".join([ch, start, end, coverage, '\n'])
        outfile.write(line)

