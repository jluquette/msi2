#!/usr/bin/env python

"""
Several computations to fit models to polymerase stutter data derived from
haploid chromosomes.  Should also support optional plotting routines to
inspect fits by eye.
"""

import sys
import rpy
from math import ceil, log
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


def estimate_mean_and_variance(binned_data):
    """Given a histogram of observed data (binned or not), estimate the mean
    and variance as functions of major allele length, then use linear
    regression to model the mean and variance.  These fitted curves will be
    used to construct PMFs for each major allele length."""

    # mean() and var() transform the allele as a -> abs(a) - 1 so a single
    # unit gain is 0 and a single unit loss is also 0.
    def mean(a, f):
        N = sum(f)
        return sum(f*(abs(a) - 1.0) for a, f in zip(a, f)) / N if N > 0 else None

    def var(a, f):
        N = sum(f)
        if N > 1:
            mu = mean(a, f)
            return sum(f*(abs(a) - 1.0 - mu)**2 for a, f in zip(a, f)) / (N - 1)
        else:
            return None

    tot_ests = { 'major': [], 'mean': [], 'var': [], 'count': [] }
    gain_ests = { 'major': [], 'mean': [], 'var': [], 'count': [] }
    loss_ests = { 'major': [], 'mean': [], 'var': [], 'count': [] }
    for major_allele, data in binned_data.iteritems():
        tot_ests['major'].append(major_allele)
        tot_ests['mean'].append(mean(a=data.keys(), f=data.values()))
        tot_ests['var'].append(var(a=data.keys(), f=data.values()))
        tot_ests['count'].append(sum(data.values()))

        gains = dict((a, f) for a, f in data.items() if a > 0)
        gain_ests['major'].append(major_allele)
        gain_ests['mean'].append(mean(a=gains.keys(), f=gains.values()))
        gain_ests['var'].append(var(a=gains.keys(), f=gains.values()))
        gain_ests['count'].append(sum(gains.values()))

        losses = dict((a, f) for a, f in data.items() if a < 0)
        loss_ests['major'].append(major_allele)
        loss_ests['mean'].append(mean(a=losses.keys(), f=losses.values()))
        loss_ests['var'].append(var(a=losses.keys(), f=losses.values()))
        loss_ests['count'].append(sum(losses.values()))

    return tot_ests, gain_ests, loss_ests


def predict_lm(lm, xs, fn):
    coef, intercept = lm['coefficients'].values()
    return [ coef*fn(x) + intercept for x in xs ]


def plot_mean_and_variance(unit, event, ests, model):
    """ests - dict(major: major alleles, mean: means, var: variances
    model   - a tuple(mean model, variance model)
    TODO: add error bars"""
    x_axis = range(min(ests['major']), max(ests['major']) + 1)
    mean_model, var_model = model
    rpy.r.plot(y=[ na(x) for x in ests['mean'] ],
        x=ests['major'], xlab="Major allele", ylab="Freq",
        main="%s, %s: mean vs. major allele" % (unit, event))
    rpy.r.lines(predict_lm(mean_model, x_axis, log), x=x_axis, col="blue")

    rpy.r.plot(y=[ na(var) for var in ests['var'] ],
        x=ests['major'], xlab="Major allele", ylab="Freq",
        main="%s, %s: variance vs. major allele" % (unit, event))
    rpy.r.lines(predict_lm(var_model, x_axis, log), x=x_axis, col="blue")


def plot_binned_data(file_name, binned_data, plots_per_row=6, interactive=False):
    vert_rows = int(ceil(len(binned_data) / float(plots_per_row)))
    with PlotContext(interactive, rpy.r.png, filename=file_name,
        width=plots_per_row*350, height=350*vert_rows):
        rpy.r.library('plotrix')
        rpy.r.layout(rpy.r.matrix(
            range(1, vert_rows*plots_per_row + 1), ncol=plots_per_row, byrow=True)
        )
        rpy.r.par(mar=[4,2,2,1])
        for major_allele, data in binned_data.items():
            rpy.r.barp(height=[ float(x)/sum(data.values()) for x in data.values() ],
                       x=[ major_allele + a for a in data ],
                       ylog=True,
                       main="Major allele %d" % major_allele)


def model_mean_and_variance(meanvar_ests):
    """Regression models of mean and var as functions of major allele len.
    NOTE: since alleles are already normalized to the major allele (e.g.,
    allele len=0 is the major allele), we're modeling error of the means
    and variances.
    LATER: Use several regression formulae to see how things look and choose
    the best fit?  For now, the log regression seems (simply by eye) to be
    the better fit."""
    # Weights are just the number of observed sites for each majro allele
    weights = meanvar_ests['count']
    adjmean = [ na(mean) for mean in meanvar_ests['mean'] ]
    adjvar = [ na(var) for var in meanvar_ests['var'] ]
    lmdata = rpy.r.data_frame(major=meanvar_ests['major'], mean=adjmean, var=adjvar)
    meanmodel = rpy.r.lm(rpy.r("mean ~ log(major)"), data=lmdata,
        weights=weights)
    varmodel = rpy.r.lm(rpy.r("var ~ log(major)"), data=lmdata,
        weights=weights)
    return meanmodel, varmodel


def model_error(model_counts, bin=2, save_rdata=None, interactive=False):
    """Given the error profile of the X and Y chromosome, derive an error
    model for stutter introduced into the library during MDA.  If save_rdata
    is a file name, convert data into R readable format and save it in the
    RData format."""

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
    binned_data = defaultdict(dict)
    unbinned_data = defaultdict(dict)
    for unit, counts_by_major_allele in model_counts.iteritems():
        for major_allele, data in counts_by_major_allele.iteritems():
            major_allele = int(major_allele)
            # bin_data with bin extension=0 is the unbinned data for that allele
            unbinned_data[unit][major_allele] = \
                bin_data(major_allele, 0, counts_by_major_allele)
            # NB. 'data' is not the right thing to pass
            binned_data[unit][major_allele] = \
                bin_data(major_allele, bin, counts_by_major_allele)
        #plot_binned_data("%s_binned_hists.png" % unit, binned_data[unit],
                         #interactive=interactive)
        rpy.r.assign('%s.raw.data' % unit, unbinned_data[unit])
        rpy.r.assign('%s.binned.data' % unit, binned_data[unit])

    # Step 2. with our binned profiles, we can determine the empirical mean
    # and variance for errors for each major allele size.  This needs to be
    # done separately for gains and losses since the final piecewise PMF will
    # be constructed with two different negbin models--one for gain and one
    # for loss.
    model = defaultdict(dict)
    for unit, bdata in binned_data.iteritems():
        tot_ests, gain_ests, loss_ests = \
            estimate_mean_and_variance(bdata)
        rpy.r.assign("%s.tot.ests" % unit, tot_ests)
        rpy.r.assign("%s.gain.ests" % unit, gain_ests)
        rpy.r.assign("%s.loss.ests" % unit, loss_ests)

        print("unit=" + unit)
        model[unit]['gain'] = model_mean_and_variance(gain_ests)
        model[unit]['loss'] = model_mean_and_variance(loss_ests)
        rpy.r.assign("%s.gain.model" % unit, model[unit]['gain'])
        rpy.r.assign("%s.loss.model" % unit, model[unit]['loss'])

        with PlotContext(interactive, rpy.r.png,
            filename="%s_mean_and_var.png" % unit, width=1000, height=1000):
            rpy.r.layout(rpy.r.matrix([1,2,3,4], ncol=2))
            rpy.r.par(mar=[3, 2, 2, 1])
            plot_mean_and_variance(unit, "gain", gain_ests, model[unit]['gain'])
            plot_mean_and_variance(unit, "loss", loss_ests, model[unit]['loss'])
 
    if save_rdata:
        rpy.r.save_image(file=save_rdata)


def fit(model_counts, bin, save_rdata=None, interactive=False):
    """Fit models to the data and return a dictionary of PMFs parameterized
    by repeat type (mono, di, tri, tetra) and major allele length.  These PMFs
    aim to describe the shape of a single peak."""

    model_error(model_counts, bin, save_rdata=save_rdata, interactive=interactive)
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


class PlotContext(object):
    """Request permission to create an R plot.  In interactive mode, 
    plot to the screen over X11 forwarding and wait for the user to hit the
    ENTER key to move to the next plot; in non-interactive mode, plot to file.
    Implements python's context guard for using "with" operator."""

    def __init__(self, interactive, plotfn, **kwargs):
        self.interactive = interactive
        self.plotfn = plotfn
        self.plotfn_kwargs = kwargs

    def __enter__(self):
        if not self.interactive:
            self.plotfn(**self.plotfn_kwargs)
        return self

    def __exit__(self, type, value, traceback):
        if not self.interactive:
            rpy.r.dev_off()
        else:
            print("ENTER for the next plot or Q to quit.")
            if sys.stdin.read(1).lower() == "q":
                exit(1)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('error_count_file',
        help='JSON file containing model parameters estimated by model_error.py')
    parser.add_argument('--save-rdata', default=None, metavar='FILE_PATH',
        help='Save various bits of data used to fit the model in an RData file.')
    parser.add_argument('--interactive', default=False, action='store_true',
        help='Do not save R plots, instead present them to the user and ' \
             'wait for the user to hit ENTER to proceed to the next plot.')
    parser.add_argument('--bin', default=2, metavar='INT', type=int,
        help='Bin observed data from major alleles within +/- INT units to '
             'smooth out the mean and variance profiles.')
    args = parser.parse_args()

    with open(args.error_count_file, 'r') as f:
        error_profile = simplejson.load(f)
    model_counts = error_profile['counts']
    fit(model_counts, args.bin, args.save_rdata, args.interactive)
