#!/usr/bin/env python

import sys

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('usage: %s candidate_table locus_chr locus_start')
        exit(1)

    locus_chr = sys.argv[2]
    locus_start = sys.argv[3]

    f = open(sys.argv[1], 'r')
    first_line = True
    for line in f:
        fields = map(str.strip, line.split('\t'))

        if first_line:
            colnames = fields
            first_line = False
            continue

        info = dict(zip(colnames, fields))
        if info['chr'] == locus_chr and info['start'] == locus_start:
            for x in zip(colnames, fields):
                print('%25s   %s' % (x[0],x[1]))
            # Done.
            exit(0)
