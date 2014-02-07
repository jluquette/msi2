#!/usr/bin/env python
"""Try to join candidates from Tae-min's pipeline to lobSTR's calls.

I'm not sure yet whether a missing call from lobSTR should lower our
confidence."""

import sys
from bx.intervals.intersection import IntervalTree, Interval
from collections import defaultdict


def in_tree(tree, start, end):
    """Return true if (`start`, `end`) overlaps any interval in `tree`."""
    return len(tree.find(int(start), int(end))) > 0


def build_lobstr_tree(f):
    """f is an open file handle to a lobSTR genotypes.tab file"""
    # Each sample has a dictionary of trees, one per chromosome
    colnames = []
    #chrs = defaultdict(lambda: IntervalTree())
    chrs = {}
    for line in f:
        # Skip comments
        if line.strip()[0] == '#':
            continue

        fields = map(str.strip, line.split('\t'))

        # Skip the header
        if fields[0] == 'chr':
            colnames = fields
            continue

        # When we want to add lobSTR data to the tree, use other= to add
        # the whole line.
        if fields[0] not in chrs.keys():
            chrs[fields[0]] = IntervalTree()
        chrs[fields[0]].insert(int(fields[1]), int(fields[2]), zip(colnames, fields))

    return chrs 


def in_lobstr(lobstr_trees, sample, chrom, start, end):
    """Return True if the (chrom, start, end) interval intersects with a
    lobSTR interval for the specified `sample`."""
    try:
        return in_tree(lobstr_trees[sample][chrom], start, end)
    except KeyError:
        warnings.warn("no tree for sample=%s, chrom=%s" % (sample, chrom),
            RuntimeWarning)
        return False


# For each sample, read in lobSTR's calls and build an interval tree for
# quick lookups.
def read_lobstr(samples):
    trees = {}
    for s in samples:
        with open('lobstr/%s.genotypes.tab' % s, 'r') as f:
            sys.stderr.write('building tree for ' + s + '\n')
            trees[s] = build_lobstr_tree(f)

    return trees

    
# headers should be in the table, but the format should be roughly:
#   chr start end (zygosity, support)*N
# where there are N pairs of zygosity and support columns.  Zygosity
# can be one of NA, hom, het or nocall.  Support is NA or a list of
# x:y values, where x is the repeat unit length and y is the number
# of reads supporting that repeat unit length.
colnames = []
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('usage: %s tab_separated_table')
        exit(1)

    #sys.setrecursionlimit(10000)  # For IntervalNode
    samples = [ 'heart_bulk', 'neuron_2', 'neuron_3', 'neuron_51', 'neuron_6']
    trees = read_lobstr(samples)

    # read in the results from Tae-min's sputnik pipeline
    tsv_f = open(sys.argv[1].strip(), 'r')
    first_line = True
    for line in tsv_f:
        fields = map(str.strip, line.split('\t'))
    
        # Handle the header: determine which column indexes contain the data
        # we care about.
        if first_line:
            colnames = fields + [ 'in.lobstr.' + s for s in samples ]
            print('\t'.join(colnames))
            first_line = False
            continue
    
        d = dict(zip(colnames, fields))
        found = [ str(in_lobstr(trees, s, 'chr'+d['chr'], d['start'], d['end']))
                for s in samples ]
        print('\t'.join(fields + found))
