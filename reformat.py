#!/usr/bin/env python

"""
reformat.py - Takes input from Tae-min's SPUTNIK pipeline and reformats
the read count field and calls zygosity for each locus.
"""

import sys
from collections import defaultdict

if len(sys.argv) != 2:
    print("Usage: %s input_file" % sys.argv[0])
    exit(1)

# support is a dictionary of (unit length, number of supporting reads) pairs
# min_reads_for_inclusion: minimum number of supporting reads required to
# consider a unit length in the zygosity determination.
# the simple model is: if there are 0 unit lengths with enough supporting
# reads then no call is made; else if there is only 1 unit length passing
# the inclusion criterion, then call homozygous; else heterozygous.
def call_zygosity(supp_freq, min_reads_for_inclusion=1):
    alleles = [ k for k, v in supp_freq.items() if v >= min_reads_for_inclusion ]

    if len(alleles) == 0:
        return "nocall"
    elif len(alleles) == 1:
        return "hom"
    else:
        return "het"


# line format is tab-delimited: index, chrom, start, end, read count list
# read count list is a string of numbers separated by commas.  The numbers
# are the repeat unit lengths; each number represents a single read.  So a
# string like 7,7,7,8 means 3 reads support an allele with 7 repeat units
# and a single read supports an allele with 8 repeat units.
# We drop the 'index' field.
f = open(sys.argv[1], 'r')
out_headers = [ 'chr', 'start', 'end', 'zygosity', 'support.freq' ]
for line in f:
    fields = map(str.strip, line.split('\t'))

    # header row
    if fields[0] == 'index':
        print('\t'.join(out_headers))
        continue

    if not fields[4]:
        # if the read count list is empty, there are no reads supporting
        # any repeat length allele.  skip these.
        continue

    # filter(None, seq) drops false values in seq.  Empty strings in
    # python are "falsey": http://docs.python.org/2/library/stdtypes.html#truth-value-testing
    supports = filter(None, map(str.strip, fields[4].split(',')))

    # Make a frequency dict of (repeat unit, number of supporting reads)
    supp_freq = defaultdict(int)
    for x in supports:
        supp_freq[x] += 1
    
    # throw away index (field 0)
    out_fields = fields[1:4] + [ call_zygosity(supp_freq) ] + \
        [ ' '.join('%s:%s' % (k, v) for k, v in supp_freq.items()) ]

    print('\t'.join(out_fields))
