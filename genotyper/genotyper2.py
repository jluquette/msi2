#!/usr/bin/env python

"""
Genotype loci summarized by msitools and given an error profile generated
by model_error.py.

MAJOR SHORTCOMINGS:
* This genotyper is simply not aware of subclones, which is especially 
import in bulk samples.  Extending to handle subclonality is likely a huge
project (e.g., Andrew's INTEGER paper) and I don't know how feasible it is.
* BIG PROBLEM: P_G, P_L and P_I should be parameterized by the real allele
length.  It is quite clear that if the other assumption of unit independence
is reasonable, then longer STR loci should have a higher P_G and P_L.
"""

import sys
import rpy
from math import ceil
from argparse import ArgumentParser
from strlocusiterator import STRLocusIterator
from collections import defaultdict
from itertools import product, combinations
from numpy import array, prod
from scipy.stats import binom_test
from scipy.misc import comb
from operator import itemgetter
from pprint import pprint

import simplejson

from fit_and_plot import bin_data, na

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


# Only included in itertools with python >= 2.7
def combinations_with_replacement(iterable, r):
    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    pool = tuple(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield tuple(pool[i] for i in indices)


def read_prob(obs, real, unit, params, debug=False):
    """The basic model using a geometric distribution for gains and a
    binomial distribution for losses.  Determines the probability of
    seeing a single read with STR length obs when the postulated real
    allele length is `real`."""
    if debug:
        print("read_prob(%d, %d, %d, %s)" % (obs, real, unit, str(params)))
    p_G, p_L, p_I, p_Gu, p_Lu = \
        itemgetter('p_G', 'p_L', 'p_I', 'p_Gu', 'p_Lu')(params)
    if obs > real: # This read is a unit gain.  Geometric model.
        n_units_gained = (obs - real) / unit
        geom_prob = (1 - p_Gu) * (p_Gu ** n_units_gained)
        # Total prob is the geometric prob conditioned on the gain event.
        if debug:
            print("geometric prob=" + str(p_G * geom_prob))
        return p_G * geom_prob
    elif obs < real: # This read is a unit loss.  Binomial model.
        n_units_lost = (real - obs) / unit
        binom_prob = \
            comb(real/unit, n_units_lost, exact=True) * \
            (p_Lu ** n_units_lost) *  \
            ((1 - p_Lu)**(real - n_units_lost))
        # Total prob is the binomial prob conditioned on the loss event.
        if debug:
            print("binomial prob=" + str(p_L*binom_prob))
        return p_L * binom_prob
    else: # Identity event.
        if debug:
            print("identity prob=" + str(p_I))
        return p_I

def lib_prob(obs, real, unit, pmf, debug=True):
    """The basic model using a geometric distribution for gains and a
    binomial distribution for losses.  Determines the probability of
    seeing a single read with STR length obs when the postulated real
    allele length is `real`."""
    return pmf(abs(obs - real) / unit)

def prod(x):
    """Compute the product of all values in x."""
    y = 1
    for z in x:
        y *= z
    return y

def single_peak(real_allele, alleles, unit, params):
    """Given a single peak model centered at `real_allele`, what is the
    probability of seeing the alleles with the observed support."""
    return (real_allele, real_allele), \
        prod(one_prob(a, real_allele, unit, params) ** data[0] for a, data in alleles.iteritems())


def double_peak(allele1, allele2, alleles, unit, params):
    """Given a double peak model with peaks at `real_allele1` and
    `real_allele2`, what is the probability of seeing the supplied alleles
    with the observed support reads.  Assumes that both underlying alleles
    are equally as likely to produce reads."""
    return prod((one_prob(a, real_allele1, unit, params)/2 +
                 one_prob(a, real_allele2, unit, params)/2) ** data[0]
                for a, data in alleles.iteritems())
 

def read_density(allele1, allele2, alleles, unit, pmf):
    """Compute the probability of each observed allele (ignoring #supp reads)
    given the two assumed real alleles `allele1` and `allele2`, then compute
    the probability of the read distribution over the real alleles assuming
    they should be evenly distributed."""
    allele_prob = prod(pmf(abs(a - allele1)/unit) / 2.0 +
                       pmf(abs(a - allele2)/unit) / 2.0
                       for a in alleles)

    for a in alleles:
        x = abs(a - allele1) / unit
        print("pmf(%d)=%f" % (x, pmf(x)))
        x = abs(a - allele2) / unit
        print("pmf(%d)=%f" % (x, pmf(x)))

    # I drew N=total depth reads from an infinite population that is
    # distributed with ratio (1/2)-epsilon between real alleles 1 and 2
    # and epsilon/(#alleles-2) for all other alleles.  Multinomial.
    # For version 1, just look at distribution of reads for the two real
    # alleles (or better yet, to easily consider nearby alleles, for each
    # allele, assign it to the closest in allele{1,2} and sum over each
    # group.  split reads in half if an allele is equidistant from both
    # allele1 and allele2.
    allele1_support = alleles.get(allele1, [0])[0]
    allele2_support = alleles.get(allele2, [0])[0]
    '''
    # XXX: Bad idea.  This causes alleles with next to no support to be given
    # support from nearby alleles with strong signals.  E.g.:
    #  8:1 9:1 10:18 11:1 is called 9/11 het because the 10 allele's reads are
    # split evenly between 9 and 11, causing the read probability to miss this
    # obvious nonhet.
    for a, data in alleles.iteritems():
        print("a=%s, data=%s, a1=%d, a2=%d" % (a, str(data), allele1, allele2))
        nreads = data[0]
        a1_dist = abs(a - allele1)
        a2_dist = abs(a - allele2)
        if a1_dist == a2_dist:
            print("a1_dist == a2_dist")
            allele1_support += nreads / 2.0
            allele2_support += nreads / 2.0
        elif a1_dist < a2_dist:
            print("a1_dist < a2_dist")
            allele1_support += nreads
        else:
            print("a1_dist > a2_dist")
            allele2_support += nreads
    '''
    allele1_support = int(allele1_support)
    allele2_support = int(allele2_support)
    tot_support = allele1_support + allele2_support
    print('allele1_support=' + str(allele1_support))
    print('allele2_support=' + str(allele2_support))
    # Min ensures we're measuring the probability of the most extreme event.
    # _p_binom is not symmetric in a,b if size=a+b
    # _d_binom is symmetric
    read_prob = rpy.r.dbinom(allele1_support, size=tot_support, prob=0.5)

    print('allele_prob=%f, read_prob=%f' % (allele_prob, read_prob))
    return allele_prob * read_prob
 

def gtcall(gts_with_probs, reflen):
    """Given a list of genotypes with probs. ((a1, a2), prob), return a call
    for the best genotype (highest prob) of (ref|het|hom, prob., gt)."""
    best_alleles, best_prob = \
        sorted(gts_with_probs, key=lambda x: x[1], reverse=True)[0]
    if len(gts_with_probs) > 1:
        pprint(gts_with_probs)
    if all(best_alleles[0] == a for a in best_alleles):
        call = 'ref' if best_alleles[0] == reflen else 'hom'
    else:
        call = 'het'
    return (call, best_prob, '/'.join(str(a) for a in best_alleles))


def genotype_locus(alleles, reflen, unit, pmf, model=read_density):
    """Consider all possible diploid genotypes, determine probabilities for
    each, and sort by highest prob."""
    #params = model_params[unit]
    unit = unit_to_int[unit]

    all_diploid_gts = combinations_with_replacement(alleles.keys(), 2)
    gts_with_probs = \
        [ ((allele1, allele2), model(allele1, allele2, alleles, unit, pmf))
          for allele1, allele2 in all_diploid_gts ]

    return gtcall(gts_with_probs, reflen)

def fit_negbin(hist):
    """Return a PMF function that will compute probability of a given allele."""
    # View negative alleles (losses) as positive values because negbin is
    # only defined on positive values.  Shift all gains and losses by 1 toward
    # the center so that mean values are closer to 0, so that the negbin
    # property variance > mean holds.
    norm_hist = dict((abs(allele) - 1, freq) for allele, freq in hist.items())
    #product = [ allele*freq for allele, freq in norm_hist.items() ]
    # Get the mean and variance
    mu = rpy.r.mean(rpy.r.rep(norm_hist.keys(), norm_hist.values()))
    var = rpy.r.var(rpy.r.rep(norm_hist.keys(), norm_hist.values()))
    # The dispersion parameter to negbin is: var = mu + mu^2 * (1/dispersion)
    #   => dispersion = mu^2 / (var - mu)
    print("got hist: " + str(hist))
    print("mu=%f, var=%f" % (mu, var))
    # PMF doesn't exist.  Loci with no matching PMF will be nocalls in the
    # final output.  Maps probability to 0 for any length.
    if var <= mu:
        return lambda x: 0
    disp = mu**2 / (var - mu)
    return lambda x: rpy.r.dnbinom(abs(x) - 1, mu=mu, size=disp)
    
    
def fit_piecewise_negbin(hist):
    """Return a PMF constructed with a piecewise negative binomial model."""
    total_count = sum(hist.values())
    gains = dict((allele, freq) for allele, freq in hist.items() if allele > 0)
    losses = dict((allele, freq) for allele, freq in hist.items() if allele < 0)
    p_id = float(hist[0]) / total_count
    p_gain = float(sum(gains.values())) / total_count
    p_loss = float(sum(losses.values())) / total_count
    gain_negbin = fit_negbin(gains)
    loss_negbin = fit_negbin(losses)

    def piecewise_pmf(x):
        if x == 0:
            return p_id
        elif x > 0:
            return p_gain * gain_negbin(x)
        else:
            return p_loss * loss_negbin(x)

    return piecewise_pmf


def plot_piecewise_pmf(hist, major_allele, pmf, xs=None):
    """Plot the derived PMF vs. the histogram of observed data."""
    if xs is None:
        xs = sorted(hist.keys())

    hist_counts = sum(hist.values())
    obs_y = [ float(hist[x])/hist_counts for x in hist.keys() ]
    fit_y = [ pmf(x) for x in xs ]

    rpy.r.barp(obs_y, x=hist.keys(), xlim=[min(xs), max(xs)], ylog=True,
        main="Piecewise PMF for major allele=%d" % major_allele,
        xlab="Minor allele (units)",
        ylab="Empirical frequency vs. model fit")
    rpy.r.points(fit_y, x=xs, type='b')


def model_mean_and_variance(binned_data):
    """Given a histogram of observed data (binned or not), estimate the mean
    and variance as functions of major allele length, then use linear
    regression to model the mean and variance.  These fitted curves will be
    used to construct PMFs for each major allele length."""

    emp_values = []  # tuple(major_allele, mean, var, gmean, gvar, lmean, lvar)
    for major_allele, data in binned_data.iteritems():
        N_tot = sum(freq for _, freq in data.iteritems())
        mean_tot = sum(freq*float(allele) for allele, freq in data.items()) / N_tot
        var_tot = sum(freq*(float(allele) - mean_tot)**2 for allele, freq in data.items()) / (N_tot - 1)

        gains = [ (allele, freq) for allele, freq in data.items() if allele > 0 ]
        N_gain = sum(freq for _, freq in gains)
        if N_gain > 1:
            mean_gain = sum(freq*float(allele) for allele, freq in gains) / N_gain
            var_gain = sum(freq*(float(allele) - mean_gain)**2 for allele, freq in gains) / (N_gain - 1)
        else:
            mean_gain = None
            var_gain = None

        losses = [ (allele, freq) for allele, freq in data.items() if allele < 0 ]
        N_loss = sum(freq for _, freq in losses)
        if N_loss > 1:
            mean_loss = sum(freq*float(allele) for allele, freq in losses) / N_loss
            var_loss = sum(freq*(float(allele) - mean_loss)**2 for allele, freq in losses) / (N_loss - 1)
        else:
            mean_loss = None
            var_loss = None

        emp_values.append((int(major_allele), mean_tot, var_tot, mean_gain, var_gain, mean_loss, var_loss))

    # These are just for getting plot boundaries
    all_means = sum([ [ tup[1], tup[3], tup[5] ] for tup in emp_values ], [])
    all_means = [ x for x in all_means if x ]  # drop Nones
    all_vars = sum([ [ tup[2], tup[4], tup[6] ] for tup in emp_values ], [])
    all_vars = [ x for x in all_vars if x ] # drop Nones

    rpy.r.layout(rpy.r.matrix([1,2], ncol=1))
    # Plots for visual inspection, move elsewhere
    ylim_mean = [ min(all_means), max(all_means) ]
    x_axis = [ tup[0] for tup in emp_values ]
    # Total mean
    rpy.r.plot([ na(tup[1]) for tup in emp_values ],
               x=x_axis, ylim=ylim_mean,
               main="Empirical mean error vs. major allele length",
               xlab="Major allele length (units)",
               ylab="Empirical mean error")
    # Mean for gain alleles
    rpy.r.lines([ na(tup[3]) for tup in emp_values ], x=x_axis, col="green")
    # Mean for loss alleles
    rpy.r.lines([ na(tup[5]) for tup in emp_values ], x=x_axis, col="red")
    rpy.r.legend("bottomleft", horiz=True, legend=["Total", "Gain", "Loss"],
        fill=["black", "green", "red"])

    # Total variance
    ylim_var = [ min(all_vars), max(all_vars) ]
    rpy.r.plot([ na(tup[2]) for tup in emp_values ],
               x=x_axis, ylim=ylim_var, log="y",
               main="Empirical variance of error vs. major allele length",
               xlab="Major allele length (units)",
               ylab="Empirical variance of error")
    # Gain variance
    rpy.r.lines([ na(tup[4]) for tup in emp_values ], x=x_axis, col="green")
    # Loss variance
    rpy.r.lines([ na(tup[6]) for tup in emp_values ], x=x_axis, col="red")
    rpy.r.legend("topleft", horiz=True, legend=["Total", "Gain", "Loss"],
        fill=["black", "green", "red"])

    # Linear regression on mean.  What sort of function should be used?
    # CAREFUL: might expect to use linear function since mean of each distn
    # should be roughly the major allele len, but these data are all
    # normalized to the major allele already--that is, 0=major allele len
    # and +/-(1,2,..) refer to deviations from the major allele.  So we're
    # modeling the error, not the mean.


def model_error(model_counts, bin=1):
    """Given the error profile of the X and Y chromosome, derive an error
    model for stutter introduced into the library during MDA."""

    # Step 1. since many major alleles do not have enough observations in
    # the data to make a good model, combine histograms from nearby major
    # alleles to increase the amount of data available while also retaining
    # some of the unique properties of the major allele.
    # Combine histograms for:
    #    | major_allele - bin | ... | major_allele + bin |
    # IMPORTANT NOTE: the bin parameter is in repeat units.  Therefore, if
    # we observe a major allele at 39 repeat units but NO alleles at 38 or
    # 40, for bin=1 the 39 histogram will get no extra data.  That is, there
    # is a distinction between the n closest alleles and alleles within n
    # repeat units.
    binned_model = defaultdict(dict)
    for unit, counts_by_major_allele in model_counts.iteritems():
        for major_allele, data in counts_by_major_allele.iteritems():
            major_allele = int(major_allele)
            binned_model[unit][major_allele] = \
                bin_data(major_allele, bin, counts_by_major_allele)
    for unit, binned_data_by_major_allele in binned_model.iteritems():
        vert_rows = ceil(len(binned_data_by_major_allele) / 8.0)
        #rpy.r.png('%s_piecewise_pmfs.png' % unit, width=8*250,
                  #height=300*vert_rows)
        #rpy.r.library('plotrix')
        #rpy.r.layout(rpy.r.matrix(range(1, vert_rows*8 + 1), ncol=8, byrow=True))
        #rpy.r.par(mar=[3,2,2,1])
        print("opening layout for %s, %d fits" % (unit, len(binned_data_by_major_allele)))
        for major_allele, binned_data in binned_data_by_major_allele.iteritems():
            #rpy.r.barp(binned_data.values(), x=[ major_allele + a for a in binned_data ])
            rpy.r.assign("%s.%d.freq" % (unit, major_allele), binned_data.values())
            rpy.r.assign("%s.%d.allele" % (unit, major_allele),[major_allele+a for a in binned_data ])
        #rpy.r.dev_off()
    rpy.r.save_image(file="hists.RData")
    exit(1)

    # Step 2. with our binned profiles, we can determine the empirical mean
    # and variance for errors for each major allele size.  This needs to be
    # done separately for gains and losses since the final piecewise PMF will
    # be constructed with two different negbin models--one for gain and one
    # for loss.
    for unit, binned_data in binned_model.iteritems():
        model_mean_and_variance(binned_data)
        sys.stdin.read(1)


def genotype(lw_params, filter_metrics_file, model_counts):
    """Genotype all loci in the STR locus file specified by lw_params.  Write
    the results to standard output and optionally write the filter metrics to
    `filter_metrics_file`."""

    model_error(model_counts)
    exit(1)

    pmfs = defaultdict(dict)
    for unit, binned_hist_by_major_allele in binned_model.iteritems():
        for major_allele, hist in binned_hist_by_major_allele.iteritems():
            print("building PMF for unit=%s, major allele=%d" % \
                (unit, major_allele))
            pmfs[unit][major_allele] = fit_piecewise_negbin(hist)

    # Plot the PMFs vs. the data used for fitting
    rpy.r.library('plotrix')
    for unit, pmfs_by_major_allele in pmfs.items():
        vert_rows = ceil(len(pmfs_by_major_allele) / 5.0)
        rpy.r.png('%s_piecewise_pmfs.png' % unit, width=1250,
                  height=300*vert_rows)
        rpy.r.layout(rpy.r.matrix(range(1, vert_rows*5 + 1), ncol=5, byrow=True))
        print("opening layout for %s, %d fits" % (unit, len(pmfs_by_major_allele)))
        for major_allele, pmf in pmfs_by_major_allele.iteritems():
            plot_piecewise_pmf(binned_model[unit][major_allele], major_allele, pmf)
        rpy.r.dev_off()
    exit(1)

    with STRLocusIterator(**lw_params) as locus_f:
        print("chr\tstart\tend\tref_len\tunit\tregion\tflank1\tflank2\tsequence\traw_alleles\tcall\tgenotype\tpval\tallele_summaries")
        for (chrom, start, end, reflen, unit, region, flank1, flank2, seq, reads) in locus_f:
            summaries = summarize_alleles(reads, 0)  # , reflen)
            alleles = " ".join("%d:%d,%.2f,%.2f" % ((k,) + v)
                               for k, v in summaries.items())
            call, pval, gt = genotype_locus(summaries, reflen, unit, pmfs[unit])
            print("%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%.5g\t%s" %
                  (chrom, start, end, reflen, unit, region, flank1, flank2,
                   seq, len(summaries), call, gt, pval, alleles))

        if filter_metrics_file:
            with open(filter_metrics_file, 'w') as f:
                for (description, value) in locus_f.filter_metrics():
                    f.write("%s\t%d\n" % (description, value))

                for (description, hist) in locus_f.hist_metrics():
                    f.write(description + "\n")
                    for k in sorted(hist.keys()):
                        f.write("%s\t%d\n" % (k, hist[k]))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('model_param_file',
        help='JSON file containing model parameters estimated by model_error.py')
    parser.add_argument('--filter-metrics-file', metavar='file', type=str,
        help='File to store metrics related to locus and read filtering.')
    parser.add_argument('--single-cell', action='store_true', default=False,
        help="Library was generated from a single cell.  Disables the " \
             "binomial model for >1 primary alleles.")
    STRLocusIterator.add_parser_args(parser)
    args = parser.parse_args()

    # Many of the command line args are STRLocusWalker parameters
    lw_params = dict(vars(args))
    del(lw_params['model_param_file'])
    del(lw_params['single_cell'])
    del(lw_params['filter_metrics_file'])

    with open(args.model_param_file, 'r') as f:
        error_profile = simplejson.load(f)
    model_counts = error_profile['counts']
    genotype(lw_params, args.filter_metrics_file, model_counts)
