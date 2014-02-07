#!/usr/bin/env python

"""
Given a SAM file and a coordinate range, filter the SAM down to only those
reads that span the ENTIRE coordinate range.

E.g.,

filter_bam.py my.sam chr10:93006702-93006778

will only produce reads that start before 93006702 AND end after 93006778.
"""

import sys

if len(sys.argv) != 3:
    print('usage: %s bam_file chr:start-stop' % sys.argv[0])
    exit(1)

try:
    chrom = sys.argv[2].strip().split(':')[0]
    start, stop = map(int, sys.argv[2].strip().split(':')[1].split('-'))
except IndexError:
    raise RuntimeError('specify coordinates as chrom:start-stop')

f = open(sys.argv[1], 'r')
for line in f:
    # Print header lines
    if line.startswith('@'):
        sys.stdout.write(line)
        continue

    fields = map(str.strip, line.split('\t'))

    this_chrom = fields[2]
    this_start = int(fields[3])
    this_len = len(fields[9])
    this_stop = this_start + this_len - 1

    if chrom == this_chrom and this_start <= start and this_stop >= stop:
        sys.stdout.write(line)
