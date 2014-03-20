#!/bin/env python2
"""
Regex benchmark
2014/03/20

A quick test to compare the performance of a couple different approaches to
pulling out spliced-leader containing reads.
"""
import re
import os
import glob
import time

def main():
    """Main"""
    # L. major test sample
    base_dir = "/cbcb/lab/nelsayed/raw_data/lminfectome"
    hpgl_id = "HPGL0075" # procyclic
    reads = glob.glob(os.path.join(base_dir, hpgl_id, 'processed/*.fastq'))

    # L. major SL sequence
    sl = "AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG"

    # test 1
    print("Test 1: Starting...")
    t0 = time.time()
    test1(reads, sl)
    t1 = time.time()
    print("Test 1: %0.2f minutes\n" % ((t1 - t0) / 60))

    # test 2
    print("Test 2: Starting...")
    t0 = time.time()
    test2(reads, sl)
    t1 = time.time()
    print("Test 2: %0.2f minutes\n" % ((t1 - t0) / 60))

    # test 3
    print("Test 3: Starting...")
    t0 = time.time()
    test3(reads, sl)
    t1 = time.time()
    print("Test 3: %0.2f minutes\n" % ((t1 - t0) / 60))

def test1(reads, sl, max_reads=float('inf')):
    """Test 1: single iteration regex approach"""
    regex = re.compile('|'.join(["^" + sl[-x:] for x in range(10, len(sl) + 1)]))

    counter = 0

    for filepath in reads:
        for i, entry in enumerate(readfq(open(filepath)), 1):
            # stop once we have reached desired number of reads
            if i > max_reads:
                break

            # get read sequence and id
            read_id = entry[0]
            read = entry[1]

            # check for sequence pattern in read
            match = re.search(regex, read)

            if match is not None:
                counter = counter + 1

    print("Test 1: %d matches" % counter)


def test2(reads, sl, max_reads=float('inf')):
    """Test 2: two-stage filtering (regex/regex)"""
    regex_step1 = re.compile("%s" % sl[-10:])
    regex_step2 = re.compile('|'.join(["^" + sl[-x:] for x in range(10, len(sl) + 1)]))

    counter = 0

    for filepath in reads:
        for i, entry in enumerate(readfq(open(filepath)), 1):
            # stop once we have reached desired number of reads
            if i > max_reads:
                break

            # get read sequence and id
            read_id = entry[0]
            read = entry[1]

            # check for sequence pattern in read
            match1 = re.search(regex_step1, read)

            if match1 is not None:
                match2 = re.search(regex_step2, read)

                if match2 is not None:
                    counter = counter + 1

    print("Test 2: %d matches" % counter)

def test3(reads, sl, max_reads=float('inf')):
    """Test 3: two-stage filtering (substring/regex)"""
    regex_step2 = re.compile('|'.join(["^" + sl[-x:] for x in range(10, len(sl) + 1)]))

    counter = 0

    for filepath in reads:
        for i, entry in enumerate(readfq(open(filepath)), 1):
            # stop once we have reached desired number of reads
            if i > max_reads:
                break

            # get read sequence and id
            read_id = entry[0]
            read = entry[1]

            #if match1 is not None:
            if sl[-10:] in read:
                match2 = re.search(regex_step2, read)

                if match2 is not None:
                    counter = counter + 1

    print("Test 3: %d matches" % counter)

# FASTQ parser
# source: https://github.com/lh3/readfq
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

if __name__ == "__main__":
    main()
