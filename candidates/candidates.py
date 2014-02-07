#!/usr/bin/env python

"""
select interesting allele configurations:

We have 7 samples:
  * 2 bulk DNA from heart and cortex, PCR amplified
  * 1 100 neuron batch, MDA amplified
  * 4 single neurons, MDA amplified

"Interesting" allele configurations satisfy the following criteria:
  1. The 4 single neurons are either segregated by the allele or the
     entire set of single neurons is segregated from the heart tissue.
  2. The segregating allele appears EITHER:
     a. in heart bulk tissue alone, with no evidence in cortex bulk or
        the 100 neuron batch,
     b. in BOTH the cortex bulk and 100 neuron batch but not in the
        heart bulk tissue.
     By requiring the allele to be either present or absent in BOTH
     PCR and MDA amplified cortex samples, we hope to control for
     artifacts specific to MDA.
"""

import sys
from operator import itemgetter

samples = [
    '1465-heartbulk',
    '1465-cortexbulk',
    '1465-caudate-100neurons',
    '1465-cortex-neuron2',
    '1465-cortex-neuron3',
    '1465-cortex-neuron51',
    '1465-cortex-neuron6'
]
    

if len(sys.argv) != 2:
    print('usage: %s candidate_file' % sys.argv[0])
    exit(1)

# Format is: allele1:depth,fracforward,meanmapq ... alleleN:... ...
# Return the set of alleles for which there is *any* evidence
#def parse_alleles(summary):
    #return set([ x.split(':')[0] for x in summary.split(' ') ])

import re
regex1 = re.compile('[-0-9]+:[^,]+')
def parse_allele_and_depth(summary):
    return [ regex1.match(x).group(0) for x in summary.split(' ') ]

regex2 = re.compile('([-0-9]+):[^,]+')
def parse_alleles(summary):
    return sorted(regex2.match(x).group(1) for x in summary.split(' '))

def parse_genotype(gt):
    return gt.split('/')

# Format is:
#   chrom, start, end, ref_len, unit, region
# followed by a variable number of 4-column groups:
#   #raw_alleles, call, genotype, allele_summary
f = open(sys.argv[1], 'r')
header_line = f.readline()
sys.stdout.write(header_line)
headers = [ h.strip() for h in header_line.split('\t') ]

# Columns containing genotype strings (allele_len/allele_len)
heart_gt_idx = headers.index('genotype.' + samples[0])
cortex_gt_idx = headers.index('genotype.' + samples[1])
caudate_gt_idx = headers.index('genotype.' + samples[2])
neuron_gt_idxs = [ headers.index('genotype.' + x) for x in samples[3:] ]

# Columns containing allele info (including #supporting reads)
heart_evidence_idx = headers.index('allele_summaries.' + samples[0])
cortex_evidence_idx = headers.index('allele_summaries.' + samples[1])
caudate_evidence_idx = headers.index('allele_summaries.' + samples[2])
neuron_evidence_idxs = \
    [ headers.index('allele_summaries.' + x) for x in samples[3:] ]


# Parse the file
nlines = 0
for line in f:
    nlines += 1
    if nlines % 10000 == 0:
        print('processed %d lines..' % nlines)

    fields = [ x.strip() for x in line.split('\t') ]

    # Skipping mono since they seem too unlikely
    #if fields[4] == 'mono':
        #continue

    #heart_gt = parse_genotype(fields[heart_gt_idx])
    #cortex_gt = parse_genotype(fields[cortex_gt_idx])
    #neurons100_gt = parse_genotype(fields[neurons100_gt_idx])

    # For the single cell samples, only the alleles in the genotype should
    # be considered.  Alleles not included in the final genotype are likely
    # to be errors of sequence, amplification, etc.
    neuron_gts = [ parse_genotype(fields[idx]) for idx in neuron_gt_idxs ]

    # To improve confidence in the alleles, ensure that they are witnessed
    # in the bulk cortex PCR-amplified sample.  They do NOT have to be the
    # called genotype in the bulk cortex sample--due to the fact that there
    # are many cells considered, we should not be surprised to see more than
    # 2 distinct alleles at any bulk locus.
    cortex_alleles = parse_alleles(fields[cortex_evidence_idx])

    # Only keep genotypes for which both alleles were seen in cortex bulk
    neuron_gts_in_cortex = \
        [ gt for gt in neuron_gts
          if gt[0] in cortex_alleles and gt[1] in cortex_alleles ]

    # Get rid of genotypes for which both alleles are seen in heart (not
    # necessarily CALLED in heart).
    heart_alleles = parse_alleles(fields[heart_evidence_idx])
    neuron_gts_not_in_heart = \
        [ gt for gt in neuron_gts_in_cortex
          if not (gt[0] in heart_alleles and gt[1] in heart_alleles) ]
        

    # A simple heuristic to handle allelic imbalance/dropout: for a gt to be
    # supported, make sure at least two single cells have it.  Justification:
    # If we see a segregating allele like: 0/0, 0/1, 0/1, 0/1, it's very
    # likely that the 0/0 allele was caused either by complete allele D.O. or
    # the call was influenced by heavy allele imbalance (e.g., 0:20, 1:3).
    # This should really be addressed by a multisample genotyper, but let's
    # see if this stop-gap solution can generate some interesting candidates.
    # XXX: relax this for now
    final_neuron_gts = [ '%s/%s' % tuple(gt) for gt in neuron_gts_not_in_heart
                         if neuron_gts_not_in_heart.count(gt) > 0 ]

    # Do the single cell alleles segregate the neurons?  Simple question:
    # is there more than one distinct genotype?
    if len(set(final_neuron_gts)) > 1:
        print('heart: %s\t%s' % (fields[heart_gt_idx], fields[heart_evidence_idx]))
        print('cortex: %s\t%s' % (fields[cortex_gt_idx], fields[cortex_evidence_idx]))
        print('caudate: %s\t%s' % (fields[caudate_gt_idx], fields[caudate_evidence_idx]))
        for i in range(3, len(samples)):
            print('%s: %s\t%s' % (samples[i], fields[neuron_gt_idxs[i-3]], fields[neuron_evidence_idxs[i-3]]))
        print('neuron_gts: ' + str(neuron_gts))
        print('cortex_alleles: ' + str(cortex_alleles))
        print('heart_alleles: ' + str(heart_alleles))
        print('neuron_gts_not_in_heart: ' + str(neuron_gts_not_in_heart))
        print('final_neuron_gts: ' + str(final_neuron_gts))
        print(line.strip())
        print(''.join(['-'] * 80))
