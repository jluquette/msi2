#!/usr/bin/env python

"""Assumption: the vast majority of non-ref repeat lengths in single cell
experiments are artifacts of polymerase stutter.  This script builds an
array of distances from the reference repeat length which can then be
viewed as a histogram in R.  The "reference" repeat length is estimated
as the repeat length with the most supporting reads."""

import sys
import numpy
from operator import iadd
from collections import defaultdict

def get_dists(info):
    def supp(sample):
        """Map each X:Y support list into a list of (X,Y) tuples.  X is the
        repeat unit length and Y is the number of supporting reads."""
        val = info['support.freq.' + sample]
        if val is None:
            return []
        return [ tuple(map(int, p.split(':'))) for p in val.split(' ') ]

    def dist(supps):
        """Compute distance to the reference allele; ref is the allele with
        the most support."""
        if not supps:
            return []
        max_supp = max(x[1] for x in supps)
        ref = [ x[0] for x in supps if x[1] == max_supp ][0]

        # Compute distance to ref allele; repeat each distance according
        # to the number of supporting reads for that distance.
        return reduce(iadd, [ [x[0] - ref]*x[1] for x in supps ], [])

    return reduce(iadd, map(dist, map(supp, samples)), [])




if len(sys.argv) != 2:
    print('usage: %s tab_separated_table')
    exit(1)

tsv_f = open(sys.argv[1].strip(), 'r')

# headers should be in the table, but the format should be roughly:
#   chr start end (zygosity, support)*N
# where there are N pairs of zygosity and support columns.  Zygosity
# can be one of NA, hom, het or nocall.  Support is NA or a list of
# x:y values, where x is the repeat unit length and y is the number
# of reads supporting that repeat unit length.
colnames = numpy.array([])
relevant_col_idxs = numpy.array([], dtype=numpy.int)
#samples = [ 'neuron_2', 'neuron_3', 'neuron_51', 'neuron_6']
samples = [ 'heart_bulk' ] #, 'cortex_bulk' ]
histogram = defaultdict(int)
for line in tsv_f:
    fields = numpy.array(map(str.strip, line.split('\t')))

    # Handle the header: determine which column indexes contain the data
    # we care about.
    if fields[0] == 'chr':
        for idx in range(0, len(fields)):
            for s in samples:
                if fields[idx].find(s) != -1:
                    relevant_col_idxs = numpy.append(relevant_col_idxs, [idx])
                    
        colnames = fields[relevant_col_idxs]
        continue

    
    # If we get here, the header should've already been parsed and we should
    # have a nonempty relevant_col_idxs.  Grab only the relevant fields and
    # translate any NA values to None.
    fields = map(lambda x: x if x != 'NA' else None, fields[relevant_col_idxs])
    info = dict(zip(colnames, fields))

    for d in get_dists(info):
        histogram[d] += 1


for k, v in histogram.items():
    print(str(k) + '\t' + str(v))
