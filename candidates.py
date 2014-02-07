#!/usr/bin/env python

import sys
import numpy

if len(sys.argv) != 2:
    print('usage: %s tab_separated_table')
    exit(1)

tsv_f = open(sys.argv[1].strip(), 'r')



# Return True or False if the supplied dict shows a difference in repeat
# unit length that is consistent with a lineage
# Tuples in this function are (repeat unit length, number of supporting reads)
def is_candidate(info, n_single=2, min_supp_for_allele=5):
    # Some utility functions
    def pair(sample):
        return (info['zygosity.' + sample], info['support.freq.' + sample])
    def is_hom(x): return x[0] == 'hom'
    def is_het(x): return x[0] == 'het'
    def is_supported_het(x, supp=min_supp_for_allele):
        """Determine number of repeat unit lengths that are supported by at
        least `supp` reads.  If more than 2 repeat unit lengths are supported,
        then consider `x` a supported heterozygous call."""
        if x[0] is None:
            return False

        freqs = map(tuple, [ f.split(':') for f in x[1].split(' ') ])
        supported_freqs = filter(lambda f: int(f[1]) >= supp, freqs)
        return len(supported_freqs) >= 2

    # Zip up zygosity and support calls into a tuple for easier processing
    heart_bulk = pair('heart_bulk')
    neurons = map(pair, samples[1:])

    # Look for sites that are homozygous in bulk heart DNA but heterozygous
    # in at least `n_single` single neuron samples.
    # Also filter het sites by support frequencies.  E.g., we may see:
    #     heart_bulk 8:14  neuronA 8:30 9:1  neuronB 8:22 7:2
    # This meets our hom and het criterion, but there are very few
    # supporting reads for the alternate alleles.
    # Some things to consider:
    #   1. alternate alleles that are very close in repeat unit length
    #      are likely false positives due to polymerase stutter
    #   2. alternate alleles that appear in multiple independent single
    #      cell neurons are less likely to be caused by stutter.
    supp_het_neurons = filter(is_supported_het, neurons)
    return is_hom(heart_bulk) and len(supp_het_neurons) >= n_single


# headers should be in the table, but the format should be roughly:
#   chr start end (zygosity, support)*N
# where there are N pairs of zygosity and support columns.  Zygosity
# can be one of NA, hom, het or nocall.  Support is NA or a list of
# x:y values, where x is the repeat unit length and y is the number
# of reads supporting that repeat unit length.
# We always want chr, start, end.
relevant_col_idxs = numpy.array([0,1,2], dtype=numpy.int)
colnames = numpy.array([])
samples = [ 'heart_bulk', 'cortex_bulk', 'neurons_100batch', 'neuron_2', 'neuron_3', 'neuron_51', 'neuron_6']
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
        print('\t'.join(colnames))
        continue

    
    # If we get here, the header should've already been parsed and we should
    # have a nonempty relevant_col_idxs.  Grab only the relevant fields and
    # translate any NA values to None.
    fields = map(lambda x: x if x != 'NA' else None, fields[relevant_col_idxs])
    info = dict(zip(colnames, fields))

    if is_candidate(info):
        # Translate None values back to 'NA'
        print('\t'.join(map(lambda x: 'NA' if x is None else x, fields)))
