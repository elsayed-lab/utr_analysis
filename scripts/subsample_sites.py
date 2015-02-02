#!/usr/bin/env python
"""
subsample_sites.py

Keith Hughitt
2014/12/24

Overview
--------
Subsamples reads to determine effect of depth on SL/Poly(A) acceptor site
distributions.

This can be useful, for example, when comparing samples from different species
in order to determine whether there are different rates of alternative
SL/Poly(A) usage that aren't simply a result of coverage depth differences.

Given a number of reads to sample and one or more GFFs, that number of reads
will be sampled randomly across all of the GFFs, and the distribution of
alternative site usage will be reported.

Usage
-----

python subsample_sites.py <num_reads> one.gff two.gff ... n.gff

"""
import gzip
import io
import os
import re
import sys
import pandas
import random

# input
num_reads = int(sys.argv[1])

print("Sampling %d reads..." % num_reads)

# GFFs
gffs = []
gff_cumsums = []

gff_colnames = ['chromosome', 'source', 'type', 'start', 'end', 'score',
                'strand', 'phase', 'attributes']
gff_coltypes = {'start': int, 'end':int, 'score': int}

for filename in sys.argv[2:]:
    # load file
    if filename.endswith('gz'):
        fp = gzip.open(filename)
    else:
        fp = open(filename)

    print("Parsing %s" % os.path.basename(filename))

    # filter out chromosome entries
    buffer = io.StringIO()

    lines = []

    for line in fp.readlines():
        line_str = line.decode('utf-8')
        if re.search('\tchromosome\t', line_str) is None:
            lines.append(line_str)
    buffer.writelines(lines)
    buffer.seek(0)

    # parse GFF and save as a pandas dataframe
    gff = pandas.read_csv(buffer, sep='\t', comment='#', names=gff_colnames,
                          dtype=gff_coltypes)
    gffs.append(gff)

    print("number of reads: %d" % gff.score.sum())

    if len(gff_cumsums) == 0:
        gff_cumsums.append(gff.score.cumsum())
    else:
        gff_cumsums.append(gff.score.cumsum() + gff_cumsums[-1].max())

    # clean-up
    fp.close

gff_cumsum_maxes = [x.max() for x in gff_cumsums]

# Determine distribution of acceptor sites for a sample of the size specified
genes = {}

read_source = []

for read_num in random.sample(range(max(gff_cumsum_maxes)), num_reads):
    # iterate over gffs
    for i, highest_readnum in enumerate(gff_cumsum_maxes):
        if read_num <= highest_readnum:
            break

    gff = gffs[i]
    gff_cumsum = gff_cumsums[i]

    read_source.append(i)

    # acceptor site
    site = gff[gff_cumsum <= read_num].tail(1)

    # edge case (left-most part of cdf)
    if len(site) == 0:
        site = gff.head(1)

    site_index = site.index[0]
    gene_name = re.search('Name=(.*);', 
                           site.loc[site_index].attributes).groups()[0]

    if gene_name not in genes:
        genes[gene_name] = [int(site.start)]
    else:
        if int(site.start) not in genes[gene_name]:
            genes[gene_name].append(int(site.start))

# distribution of number of acceptor sites
sites_per_gene = pandas.Series([len(v) for k,v in genes.items()])
print("================================")
print("Number of acceptor sites / gene:")
print("================================")
print(sites_per_gene.value_counts() / len(sites_per_gene) * 100)

print("================================")
print("Read sources:")
print("================================")
for x in set(read_source):
    print("%s: %d" % (os.path.basename(sys.argv[2:][x]), read_source.count(x)))

