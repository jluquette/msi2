#!/usr/bin/env python

"""Rank candidates by taking the max distance from the allele in bulk
heart DNA and all of the alleles in the single cells."""

import sys
import numpy
from operator import iadd

if len(sys.argv) != 2:
    print('usage: %s tab_separated_table')
    exit(1)

tsv_f = open(sys.argv[1].strip(), 'r')


def max_dist(info, min_supp_for_alleles=5):
    """Get the max difference between the repeat unit length in heart bulk
    and all of the repeat unit lengths in the single cell neurons."""
    ref = int(info['support.freq.heart_bulk'].split(':')[0])

    def supported_allele_lens(s):
        if s is None:
            return []
        pairs = [ p.split(':') for p in s.split(' ') ]
        return [ int(p[0]) for p in pairs if int(p[1]) >= min_supp_for_alleles ]

    neurons = [ info['support.freq.' + s] for s in samples[1:] ]
    allele_lens = reduce(iadd, map(supported_allele_lens, neurons))
    return max(abs(l - ref) for l in allele_lens)

 
# headers should be in the table, but the format should be roughly:
#   chr start end (zygosity, support)*N
# where there are N pairs of zygosity and support columns.  Zygosity
# can be one of NA, hom, het or nocall.  Support is NA or a list of
# x:y values, where x is the repeat unit length and y is the number
# of reads supporting that repeat unit length.
colnames = numpy.array([])
all_colnames = []
relevant_col_idxs = numpy.array([], dtype=numpy.int)
samples = [ 'heart_bulk', 'neuron_2', 'neuron_3', 'neuron_51', 'neuron_6']
results = []
first_line = True
for line in tsv_f:
    fields = numpy.array(map(str.strip, line.split('\t')))

    # Handle the header: determine which column indexes contain the data
    # we care about.
    if first_line:
        all_colnames = fields
        for idx in range(0, len(fields)):
            for s in samples:
                if fields[idx].find(s) != -1:
                    relevant_col_idxs = numpy.append(relevant_col_idxs, [idx])
                    
        colnames = fields[relevant_col_idxs]
        first_line = False
        continue

    
    # If we get here, the header should've already been parsed and we should
    # have a nonempty relevant_col_idxs.  Grab only the relevant fields and
    # translate any NA values to None.
    none_fields = map(lambda x: x if x != 'NA' else None, fields[relevant_col_idxs])
    info = dict(zip(colnames, none_fields))
    mdist = max_dist(info)

    
    results += [(mdist, fields)]


results = sorted(results, key=lambda entry: entry[0], reverse=True)

print('\t'.join(['max_dist'] + all_colnames.tolist()))
for entry in results:
    print('\t'.join([str(entry[0])] + entry[1].tolist()))
