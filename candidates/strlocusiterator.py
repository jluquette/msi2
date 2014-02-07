#!/usr/bin/env python

from collections import defaultdict

unit_to_int = { 'mono': 1, 'di': 2, 'tri': 3, 'tetra': 4 }

def build_reads(obslens, strands, mapqs, min_mapq):
    """Return a list of tuples describing each read a this locus:
            (observed repeat length, strand, mapq),
    Only include reads with mapq >= `min_mapq`
    """
    # "if x": don't try to int() cast an empty string, just drop it
    obs = ( int(x.strip()) for x in obslens.split(',') if x )
    strands = ( x.strip() for x in strands.split(',') if x )
    mapqs = ( int(x.strip()) for x in mapqs.split(',') if x )

    return [ read for read in zip(obs, strands, mapqs) if read[2] >= min_mapq ]


class STRLocusIterator():
    """Iterate over an STR locus file as produced by msitools.
    The line format is tab-delimited:
        chrom, start, end, unit, region, STR len list, strand list, mapq list.
    STR len list is a string of numbers separated by commas.  The numbers are
    the repeat unit lengths in bp; each number represents a single read.  So a
    string like 7,7,7,8 means 3 reads support an allele with 7bp of the STR and
    a single read supports an allele with 8bp of the STR.  The strand and mapq
    lists are similarly formatted.  Positions in the 3 lists all correspond to
    the same read.  So if strlens=7,7,8, strands=+,-,- and mapqs=60,20,40; then
    the first read has an STR len of 7bp, is on the '+' strand and has mapq=60.

    For each locus, return a tuple describing the locus as well as a list of
    tuples describing each read at the locus.  Loci and reads are subject to
    optional filtering criteria; any reads or loci failing the filter criteria
    will be quietly dropped or passed over by this iterator.  Counts for reads
    and loci dropped by the filtering criteria are tracked internally and can
    be printed at any time.
    """

    def __init__(self, filename, min_mapq=0, min_units=0,
                 max_ref_diff=float('+inf'),
                 min_supp_reads=0, x_only=False, y_only=False):
        """Filter options:
        min_mapq   -- only consider reads with mapq >= `min_mapq`
        min_units  -- only consider repeat loci with ref_units > `min_units`
        x_only     -- only consider repeat loci on the X chromosome
        y_only     -- only consider repeat loci on the Y chromosome
        """

        if x_only and y_only:
            raise RuntimeError('only one of x_only and y_only my be specified')

        self.f = open(filename, 'r')
        self.header = self.f.readline()

        self.min_mapq = min_mapq
        self.x_only = x_only
        self.y_only = y_only
        self.min_units = min_units
        self.min_supp_reads = min_supp_reads
        self.max_ref_diff = max_ref_diff

        self.total_loci = 0
        self.total_reads = 0
        self.pf_loci = 0
        self.pf_reads = 0
        self.loci_x_only_filter = 0
        self.reads_x_only_filter = 0
        self.loci_y_only_filter = 0
        self.reads_y_only_filter = 0
        self.loci_mapq_filter = 0
        self.reads_mapq_filter = 0
        self.loci_min_units_filter = 0
        self.reads_min_units_filter = 0
        self.loci_max_ref_diff_filter = 0
        self.reads_max_ref_diff_filter = 0
        self.loci_min_supp_filter = 0
        self.reads_min_supp_filter = 0
        self.chrom_hist = defaultdict(int)
        self.reflen_hist = defaultdict(int)
        self.reflen_diff_hist = defaultdict(int)
        self.nsupp_hist = defaultdict(int)
        self.nalleles_hist = defaultdict(int)
        self.strand_hist = defaultdict(int)
        self.mapq_hist = defaultdict(int)
        self.region_hist = defaultdict(int)
        self.unit_hist = defaultdict(int)
        self.n_units_hist = defaultdict(int)


    @staticmethod
    def add_parser_args(parser):
        """Add options understood by __init__ to `parser`."""
        parser.add_argument('filename', metavar='str_summary', type=str,
            help='STR summary file from msitools')
        parser.add_argument('--min-mapq', dest='min_mapq', metavar='N',
            default=0, type=int,
            help='Discard reads with mapping quality < N')
        parser.add_argument('--min-units', dest='min_units', metavar='N',
            default=1, type=int,
            help='Discard reference loci with < N repeat units')
        parser.add_argument('--min-supp', dest='min_supp_reads', metavar='N',
            default=0, type=int,
            help='Discard reference loci with < N supporting reads ' \
                 'after all read filters have been applied')
        parser.add_argument('--max-ref-diff', dest='max_ref_diff', metavar='N',
            default=float('+inf'), type=int,
            help='Discard reads that differ too greatly from the ' \
                 'reference STR length.  abs(observed len - ref len) ' \
                 '> N, for N in base pairs')
        parser.add_argument('--x-only', dest='x_only', action='store_true',
            default=False,
            help='Only consider reads from the X chromosome')
        parser.add_argument('--y-only', dest='y_only', action='store_true',
            default=False,
            help='Only consider reads from the Y chromosome')


    def __iter__(self):
        return self


    def next(self):
        while True:
            line = self.f.next()
            fields = [ x.strip() for x in line.split('\t') ]

            chrom = fields[0].replace('chr', '')
            start = int(fields[1])
            end = int(fields[2])
            ref_unit = fields[3]
            region = fields[4]

            # Need to tally unfiltered totals before any filtering occurs
            self.total_loci += 1
            raw_reads = fields[5].count(',') + 1 if fields[5] else 0
            self.total_reads += raw_reads

            if self.x_only and not chrom.endswith('X'):
                self.loci_x_only_filter += 1
                self.reads_x_only_filter += raw_reads
                continue

            if self.y_only and not chrom.endswith('Y'):
                self.loci_y_only_filter += 1
                self.reads_y_only_filter += raw_reads
                continue

            # reference STR length in bp
            reflen = end - start + 1  
            units = float(reflen) / unit_to_int[ref_unit]
            if units < self.min_units:
                self.loci_min_units_filter += 1
                self.reads_min_units_filter += raw_reads
                continue

            # No reads here
            if not raw_reads:
                continue

            # Build read tuples and filter by mapQ
            reads_mq = build_reads(fields[5], fields[6], fields[7], self.min_mapq)
            self.reads_mapq_filter += raw_reads - len(reads_mq)
            if not reads_mq:
                self.loci_mapq_filter += 1
                continue

            # Filter by difference between ref and obs alleles.  Some very
            # large differences cannot be supported by our method; other very
            # large differences are very low confidence.  In some cases, this
            # will throw away real differences (e.g., big deletions).
            reads = [ read for read in reads_mq
                               if abs(read[0] - reflen) < self.max_ref_diff ]
            self.reads_max_ref_diff_filter += len(reads_mq) - len(reads)
            if not reads:
                self.loci_max_ref_diff_filter += 1
                continue

            # Filter by total supporting reads at the locus
            if len(reads) < self.min_supp_reads:
                self.loci_min_supp_filter += 1
                self.reads_min_supp_filter += len(reads)
                continue

            # No more filters beyond this point.  Track locus statistics.
            self.pf_loci += 1
            self.pf_reads += len(reads)
            self.nsupp_hist[len(reads)] += 1
            self.chrom_hist[chrom] += len(reads)
            self.unit_hist[ref_unit] += 1
            self.n_units_hist[units] += 1
            self.reflen_hist[reflen] += len(reads)
            self.region_hist[region] += 1

            # Track read statistics
            for (obslen, strand, mapq) in reads:
                self.reflen_diff_hist[obslen - reflen] += 1
                self.strand_hist[strand] += 1
                self.mapq_hist[mapq] += 1

            # Number of unique observed alleles
            self.nalleles_hist[len(set(obslen for (obslen, x, y) in reads))] += 1

            return (chrom, start, end, ref_unit, region, reads)


    def filter_metrics(self):
        """Return a list of filter metrics.  (pf = passing filters)"""
        return [
            ('total_reads'              , self.total_reads),
            ('total_loci'               , self.total_loci),
            ('pf_reads'                 , self.pf_reads),
            ('pf_loci'                  , self.pf_loci),
            ('filter_x_only_reads'      , self.reads_x_only_filter),
            ('filter_x_only_loci'       , self.loci_x_only_filter),
            ('filter_y_only_reads'      , self.reads_y_only_filter),
            ('filter_y_only_loci'       , self.loci_y_only_filter),
            ('filter_mapq_reads'        , self.reads_mapq_filter),
            ('filter_mapq_loci'         , self.loci_mapq_filter),
            ('filter_min_units_reads'   , self.reads_min_units_filter),
            ('filter_min_units_loci'    , self.loci_min_units_filter),
            ('filter_min_supp_reads'    , self.reads_min_supp_filter),
            ('filter_min_supp_loci'     , self.loci_min_supp_filter),
            ('filter_max_ref_diff_reads', self.reads_max_ref_diff_filter),
            ('filter_max_ref_diff_loci' , self.loci_max_ref_diff_filter)
        ]


    def hist_metrics(self):
        """Return several distributions as (description, dict) tuples."""
        return [
            ("supp reads by strand", self.strand_hist),
            ("supp reads by chromosome", self.chrom_hist),
            ("supp reads by ref STR len", self.reflen_hist),
            ("supp reads by difference wrt ref STR len in bp", self.reflen_diff_hist),
            ("supp reads by mapping quality", self.mapq_hist),
            ("ref STR loci by region type", self.region_hist),
            ("ref STR loci by unit length", self.unit_hist),
            ("ref STR loci by number of units", self.n_units_hist),
            ("ref STR loci by supp reads", self.nsupp_hist),
            ("ref STR loci by num of alleles", self.nalleles_hist)
        ]
