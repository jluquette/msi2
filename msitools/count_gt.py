#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    print('usage: %s sputnik_results_file' % sys.argv[0])
    exit(1)

f = open(sys.argv[1], 'r')
for line in f:
    fields = map(str.strip, line.split(' '))
    seq = fields[-1]

    print(' '.join(fields) + ' ' + str(seq.count('GT')))
