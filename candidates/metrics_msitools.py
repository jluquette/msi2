#!/usr/bin/env python

"""
metrics_msitools.py - Make a frequency table of #supporting reads per locus.
Compute some other interesting stats as well.  Stats in brackets would be
interesting but require joining the original repeat database.
  * how many sites have evidence for multiple alleles
  * mapQ histogram
  * strand distribution
  * chrom distribution
  * STR len (in bp) distribution
  * difference between observed STR len and ref STR len
    not absolute value: negative values show obs < ref len
  * STR unit distribution: mono, di, tri, tetra
  * STR unit multiple distribution: (end-start) / unit
  * STR genomic location distribution: intergenic, intronic, exonic

Options for filtering reads:
  * haploid options: --only-x, --only-y, only use X or Y alignments, depending
    on option (both cannot be specified together).  For male samples, these
    should be haploid and thus deviations from expectation should reflect
    experimental errors.
  * maximum mapQ option: --mapq60, only use mapq60 (very high confidence)
    alignments
  * minimum repeat unit option: --min-unit (recommended=3, default=1?).  Only
    consider loci where the repeat unit (end-start+1) / unit size is greater
    than the specified value.  This could be useful to remove questionable
    loci, like 2~3 repeat units of tri or tetranucleotide repeats.
"""

import sys
from argparse import ArgumentParser
from strlocusiterator import STRLocusIterator

parser = ArgumentParser()
STRLocusIterator.add_parser_args(parser)
args = parser.parse_args()

locus_f = STRLocusIterator(**vars(args))
for (chrom, start, end, unit, region, reads) in locus_f:
    # Don't do anything, just accumulate metrics
    continue

for (description, value) in locus_f.filter_metrics():
    print("%s\t%d" % (description, value))

for (description, hist) in locus_f.hist_metrics():
    print(description)
    for k in sorted(hist.keys()):
        print("%s\t%d" % (k, hist[k]))
