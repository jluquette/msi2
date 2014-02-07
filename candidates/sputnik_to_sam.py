#!/usr/bin/env python
import sys

if len(sys.argv) != 2:
    print('usage: %s sputnik_results_file')
    exit(1)

f = open(sys.argv[1], 'r')
for line in f:
    fields = map(str.strip, line.split(' '))

    diff = int(fields[6]) - int(fields[3])
    print('\t'.join(fields[0:4] + ['1'] + fields[4:7] + [str(diff), fields[12], 'J'*len(fields[12])]))
