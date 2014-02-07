#!/usr/bin/env python

"""
Profile STR polymorphism error rates and genotype loci summarized by
msitools.  To determine error rate, assume
  (a) that the output comes from single cell experiments and
  (b) come from an XY-male specimen
So that we may assume that every locus on the X and Y chromosomes must be
mono-allelic and that every read/variant allele must be an error.

MAJOR SHORTCOMINGS:
* This genotyper is simply not aware of subclones, which is especially 
import in bulk samples.  Extending to handle subclonality is likely a huge
project (e.g., Andrew's INTEGER paper) and I don't know how feasible it is.
* BIG PROBLEM: P_G, P_L and P_I should be parameterized by the real allele
length.  It is quite clear that if the other assumption of unit independence
is reasonable, then longer STR loci should have a higher P_G and P_L.
"""

from argparse import ArgumentParser
from strlocusiterator import STRLocusIterator
from collections import defaultdict
from itertools import product, combinations
from numpy import array, prod
from scipy.stats import binom_test
from scipy.misc import comb
from operator import itemgetter
from pprint import pprint
from simplejson import dump


unit_to_int = { 'mono':1, 'di':2, 'tri':3, 'tetra':4 }

def summarize_alleles(reads, reflen):
    """`reads` is a list of (obs len, strand, mapq) tuples each describing
    a single read.
    `reflen` is the length of the STR in the reference genome."""
    # Map: observed len -> (num reads, num forward strand, sum mapq)
    metrics = defaultdict(lambda: array([0, 0, 0]))
    for obs, strand, mapq in reads:
        metrics[obs - reflen] += array([ 1, 1 if strand == '+' else 0, mapq ])

    return dict((k, (n, float(nforward) / n, float(summapq) / n))
                for k, (n, nforward, summapq) in metrics.iteritems())


def div_or_zero(x, y):
    """Convert to floating point, divide if possible, return 0 if undefined."""
    return float(x)/y if y != 0 else 0


class ErrorProfile():
    """Cumulative error count and keep a list of all witnessed error counts."""
    def __init__(self):
        self.loci = 0  # == len(self.ecs)
        #self.ecs = []
        # estimate_error really defines the contents of an ErrorCount object.
        # using an empty list ensures a 0-count object.
        self.cumulative_ec = ErrorCounts({}, 1)

    def witness(self, error_count):
        self.loci += 1
        self.cumulative_ec += error_count
        #self.ecs.append(error_count)

    def __str__(self):
        return "Raw data: %s\nProbs: %s" % \
            (str(self.cumulative_ec), str(self.estimate_model_params()))

    def estimate_model_params(self, ec=None):
        """Estimate several probabilities.  The error profile is:
            (#gain events, #loss events, #identity events, #gain successes,
            #gain attempts, #loss successes, #loss attempts)
            1. P(G), the probability of any gain event
            2. P(L), the probability of any loss event
            3. P(I), the probability of an identity event (i.e., no error)
            4. P_Gu, the probability of a single unit gain event, given a gain
                    has occurred.
            5. P_Lu, the probability of a single unit loss event, given a loss
                    has occurred."""
    
        if ec is None:
            ec = self.cumulative_ec
        self.p_G = div_or_zero(ec.gain_events, ec.gain_events + ec.loss_events + ec.identity_events)
        self.p_L = div_or_zero(ec.loss_events, ec.gain_events + ec.loss_events + ec.identity_events)
        self.p_I = div_or_zero(ec.identity_events, ec.gain_events + ec.loss_events + ec.identity_events)
        self.p_Gu = div_or_zero(ec.gain_successes, ec.gain_attempts)
        self.p_Lu = div_or_zero(ec.loss_successes, ec.loss_attempts)

        self.p_LibG = div_or_zero(ec.library_gains, ec.library_gains + ec.library_losses + ec.library_identities)
        self.p_LibL = div_or_zero(ec.library_losses, ec.library_gains + ec.library_losses + ec.library_identities)
        self.p_LibI = div_or_zero(ec.library_identities, ec.library_gains + ec.library_losses + ec.library_identities)
        self.p_LibGu = div_or_zero(ec.library_gain_successes, ec.library_gain_attempts)
        self.p_LibLu = div_or_zero(ec.library_loss_successes, ec.library_loss_attempts)

        return {
            'p_LibG': self.p_LibG,
            'p_LibL': self.p_LibL,
            'p_LibI': self.p_LibI,
            'p_LibGu': self.p_LibGu,
            'p_LibLu': self.p_LibLu,
            'p_G': self.p_G,
            'p_L': self.p_L,
            'p_I': self.p_I,
            'p_Gu': self.p_Gu,
            'p_Lu': self.p_Lu
        }


def profile_error_distn(lw_params):
    """Create a profile of likely erroneous read -> STR allele mappings.
    The method works by examining STR polymorphisms on hemizygous sex
    chromosomes (e.g., depends on the subject being male). At each locus:
      1. the primary allele is determined by the allele with the highest
         number of supporting reads.
      2. all other alleles are considered to be PCR errors.  The alleles
         are separated into gain and loss events and tabulated.  The
         computed tables are used to estimate parameters for a geometric
         gain model and a binomial loss model.
    NOTE: apart from choosing the primary allele, alleles are not weighted
    by the number of reads supporting them.  That is, we are not modeling
    error during sequencing, but rather error during library construction."""
    total_loci = 0
    errprofs = defaultdict(lambda: defaultdict(lambda: ErrorProfile()))
    lw_params['hemizygous_only'] = True  # force this, may want to allow walking over autosomes at a later point
    with STRLocusIterator(**lw_params) as locus_f:
        for (chrom, start, end, reflen, unit, region, flank1, flank2, seq, reads) in locus_f:
            # Using reflen=0 causes the summaries to contain the absolute
            # number of bp in the STR locus, not the deviation from ref.
            if total_loci % 10000 == 0:
                print("PROFILE ERROR: finished parsing %d loci.." % total_loci)
            summaries = summarize_alleles(reads, 0)
            ec = ErrorCounts(summaries, unit_to_int[unit])
            errprofs[unit][ec.real_len].witness(ec)
            total_loci += 1

    print("Total loci in error profile: %d" % total_loci)

    model_counts = defaultdict(dict)
    model_params = defaultdict(dict)
    for unit, prof_by_real_len in errprofs.items():
        for real_len, eprof in prof_by_real_len.items():
            #model_params[unit][real_len] = eprof.estimate_model_params()
            model_counts[unit][real_len] = eprof.cumulative_ec.counts()
    print("Supporting counts for parameters:")
    pprint(model_counts)
    #print("Inferred model parameters:")
    #pprint(model_params)

    #return model_counts, model_params
    return model_counts


class ErrorCounts(object):
    def __init__(self, summaries, unit_len):
        """Return an object summarizing several error modes present at this
        locus, where "real_len" is the approximation for the real allele len.
        This approximation is currently performed by selecting the allele with
        the most supporting reads.  Unrelated to the reference length.
        NOTE: can be called on summaries={} to return a 0 ErrorCounts object.
            1. number of gain events (obslen > real_len)
            2. number of loss events (obslen < real_len)
            3. number of identity events (obslen = real_len)
        --- In the error model, we will always condition gain and loss
        --- models ON the event that a gain or loss has already occurred,
        --- so we should never consider identity events when tabulating
        --- gain units and "gain unit attempts" (and likewise for loss).
            4. number of gained units
               (sum(obslen - real_len) when obslen > real_len)
            5. number of unit gain attempts (e.g., for any gain event G of n
               units, the #gained is n and #gain attempts is n+1, because the
               model is that there were n independent gain "successes" followed
               by a single gain "failure", which translates to a standard
               geometric distribution.)
            6. number of lost units (sum(real_len - obslen) when obslen < real_len)
            7. number of unit loss attempts (e.g., for any loss event L of n
               units at a site with real_len=N, the #loss "successes" is n
               and the #loss "attempts" is N, which corresponds to a binomial
               probability of n successes in N independent trials.  XXX: I'm
               still not entirely convinced that this model is correct, because
               we will be tabulating successes and attempts over sites with
               varying real_len.  However, the justification for this is the
               chance of PCR slippage causing a deleted unit should be
               dependent only on the unit itself (either length or exact seq),
               not the context (e.g., the size of the repeat containing the unit))."""
        # Sort by allele length.  Also convert basepair differences into unit
        # differences.
        # Summaries is a dict of tuples: abslen -> nreads, fracforward, meanmapq
        #self.hist = defaultdict(list)
        self.hist = defaultdict(int)
        self.libhist = defaultdict(int)
        self.loci = 1 if summaries else 0
        self.unit_len = unit_len

        # We don't care about anything but the number of reads now.
        alleles = \
            sorted([ (abslen / unit_len,) + data
                    for abslen, data in summaries.iteritems() ], 
                key=lambda x: x[1], reverse=True)  # x[1] = #reads
    
        # Estimate the real allele as the one with the most supp reads
        # (allow for empty alleles list)
        real_len = alleles[0][0] if alleles else 0
        self.real_len = real_len   # Keep a copy

        # Total number of reads at this locus
        tot_reads = sum(nreads for a, nreads, fracforward, meanmapq in alleles)
            
        # Now that we know the real allele, record the frequencies of
        # variation relative to the real allele.
        for a, nreads, fracforward, meanmapq in alleles:
            #self.hist[a - real_len].append(float(nreads)/tot_reads)
            self.hist[a - real_len] += nreads
            self.libhist[a - real_len] += 1
    
        '''
        gains = [ a for a in alleles if a[0] > real_len ]
        losses = [ a for a in alleles if a[0] < real_len ]
        # Should be = nreads[0][1]
        self.library_identities = 1 if alleles else 0
        self.identity_events = sum([ a[1] for a in alleles if a[0] == real_len ])
        
        # Gain model is geometric.  Each observed allele is the result of n
        # independent single unit extensions finalized by a single extension
        # failure.  This is NOT weighted by the number of reads supporting
        # each event.
        self.library_gains = len(gains)
        self.library_gain_successes = sum([ (g[0] - real_len) for g in gains ])
        self.library_gain_attempts = self.library_gains + self.library_gain_successes
        self.gain_events = sum([ g[1] for g in gains ])
        self.gain_successes = sum([ (g[0] - real_len)*g[1] for g in gains ])
        self.gain_attempts = self.gain_successes + self.gain_events
    
        # Loss model is binomial, but here we are only trying to estimate the
        # probability of unit deletion, which I postulate is dependent only on
        # repeat unit, not the size of the STR locus.  For each loss allele with
        # n deleted units, we consider the real length of the locus real_len as
        # the number of possible deletions, of which n were successful.
        self.library_losses = len(losses)
        self.library_loss_successes = sum([ (real_len - l[0]) for l in losses ])
        self.library_loss_attempts = real_len * len(losses)
        self.loss_events = sum([ l[1] for l in losses ])
        self.loss_successes = sum([ (real_len - l[0])*l[1] for l in losses ])
        self.loss_attempts = sum([ real_len*l[1] for l in losses ])
        '''


    def __add__(self, other):
        self.loci += other.loci
        '''
        self.library_identities += other.library_identities
        self.library_gains += other.library_gains
        self.library_gain_attempts += other.library_gain_attempts
        self.library_gain_successes += other.library_gain_successes
        self.library_losses += other.library_losses
        self.library_loss_attempts += other.library_loss_attempts
        self.library_loss_successes += other.library_loss_successes
        self.library_losses += other.library_losses
        self.identity_events += other.identity_events
        self.gain_events += other.gain_events
        self.gain_attempts += other.gain_attempts
        self.gain_successes += other.gain_successes
        self.loss_events += other.loss_events
        self.loss_attempts += other.loss_attempts
        self.loss_successes += other.loss_successes
        '''
        for k, v in other.hist.items():
            self.hist[k] += v
        for k, v in other.libhist.items():
            self.libhist[k] += v
        return self

    def counts(self):
        return {
            'Loci': self.loci,
            'LibHist': self.libhist,
            'Hist': self.hist
        }
        '''
            'LibI': self.library_identities,
            'LibG': self.library_gains,
            'LibGa': self.library_gain_attempts,
            'LibGs': self.library_gain_successes,
            'LibL': self.library_losses,
            'LibLa': self.library_loss_attempts,
            'LibLs': self.library_loss_successes,
            'I': self.identity_events,
            'G': self.gain_events,
            'Ga': self.gain_attempts,
            'Gs': self.gain_successes,
            'L': self.loss_events,
            'La': self.loss_attempts,
            'Ls': self.loss_successes
        }
        '''

    def __str__(self):
        return ", ".join("%s=%d" % (lab, n) for lab, n in self.counts().items())


def strunit(unit, profs):
    return "\t%s=list(\n%s\n\t)" % \
        (unit,
         ",\n".join(strprof(real_len, prof) for real_len, prof in profs.items()))

def strprof(real_len, prof):
    return "\t\t'%d'=list(\n\t\t\t%s\n\t\t)" % (real_len, strhist(prof['LibHist']))

def strhist(h):
    alleles = h.keys()
    freqs = [str(h[a]) for a in alleles ]
    return "'alleles'=c(%s),\n\t\t\t'freqs'=c(%s)" % \
        (", ".join(str(a) for a in alleles), ", ".join(freqs))

def make_r_plot_programs(f, model_counts):
    # Step 1) convert the entire counts structure into:
    # list of (mono|di|tri|tetra) ->
    #   list of real_len ->
    #       list of (allele len, freq) tuples
    f.write("histograms <- list(\n")
    f.write(",\n".join(strunit(unit, profs) for unit, profs in model_counts.items()))
    f.write("\n)")
    # Step 2) compute negative binomial fits to the data
    for unit, profs in model_counts.items():
        for real_len, eprof in profs.items():
            pass
    # Step 3) plot fitted curves against data as a simple sanity check

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('error_count_file',
        help='Save model parameters as a JSON structure in this file.')
    parser.add_argument('--filter-metrics-file', metavar='file', type=str,
        help='File to store metrics related to locus and read filtering.')
    parser.add_argument('--r-plot-program', metavar='file',
        help='Write an R program to plot stutter histograms against the ' \
             'derived model parameters.')
    STRLocusIterator.add_parser_args(parser)
    args = parser.parse_args()

    # Many of the command line args are STRLocusWalker parameters
    lw_params = dict(vars(args))
    del(lw_params['filter_metrics_file'])
    #del(lw_params['single_cell'])
    del(lw_params['error_count_file'])
    del(lw_params['r_plot_program'])

    #model_counts, model_params = profile_error_distn(lw_params)
    model_counts = profile_error_distn(lw_params)
    with open(args.error_count_file, 'w') as f:
        #dump({ 'counts': model_counts, 'params': model_params }, f)
        dump({ 'counts': model_counts }, f)
    if args.r_plot_program:
        with open(args.r_plot_program, 'w') as f:
            make_r_plot_programs(f, model_counts)
