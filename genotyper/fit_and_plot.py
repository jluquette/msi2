#!/usr/bin/env python

def bin_data(major_allele, bin, data):
    """Return a binned histogram for the supplied list of alleles."""
    binned_hist = defaultdict(int)
    for allele_to_bin in range(major_allele - bin, major_allele + bin + 1):
        # Sum up all of the hists into a single binned hist
        hist = data.get(str(allele_to_bin), {'LibHist':{}})['LibHist']
        for allele, freq in hist.items():
            binned_hist[int(allele)] += freq

    return dict(binned_hist)


def build_profiles_for_major_alleles(data, bins):
    profiles = defaultdict(list)
    # Bin value allows binning frequency histograms using major alleles in
    # a window of +/-(bin) to increase the amount of data available for each
    # histogram.
    for bin in bins:
        for major_allele in data:
            # Bin individual histograms by the +/-(bin) alleles
            major_allele = int(major_allele)
            allele_freq = bin_data(major_allele, bin, data).items()

            # Compute mean and variance when they both exist
            total = sum(freq for allele,freq in allele_freq)
            if total > 1:
                mean = sum(freq*float(allele) for allele,freq in allele_freq) / total
                var = sum(freq*(float(allele) - mean)**2 for allele,freq in allele_freq) / (total - 1)
            else:
                print("GOT total=%d FOR BIN=%d, MAJOR ALLELE=%d, UNIT=%s" % \
                    (total, bin, major_allele, unit))
                mean = None
                var = None

            profiles[bin] += [ (major_allele, mean, var) ]

        # Necessary for plotting with type=line in R
        profiles[bin] = sorted(profiles[bin], key=lambda x: x[0])

    return profiles


# Convert python None values to R NA values
def na(x): return rpy.r.NA if x is None else x

def plot_mean_and_var_against_major_allele(profiles, unit):
    # Plot mean (tup[1]) of distn against major allele len (tup[0])
    # Black for bin=0 and bin=max
    rpy.r.png("%s_mean_vs_major_allele.png" % unit, height=400, width=800)
    cols = [ 'black' ] + rpy.r.rainbow(len(profiles) - 2) + [ 'black' ]
    #rpy.r.layout(rpy.r.matrix([1,2], ncol=2))
    rpy.r.par(mar=[5, 4, 3, 1])
    for bin, col in zip(profiles, cols)[0:-1]:
        profile = profiles[bin]
        plotfn = rpy.r.plot if bin == 0 else rpy.r.lines
        plotfn(y=[ na(tup[1]) for tup in profile ], x=[ tup[0] for tup in profile ], type='l', col=col,
            xlab="Major allele length (units)",
            ylab="mean(distn of stutter)",
            main ="%s: dependence of mean on major allele length" % unit)

    # bin=len plot
    profile = profiles[max(profiles)]
    rpy.r.lines(y=[ na(tup[1]) for tup in profile ], x=[ tup[0] for tup in profile ], type='l', col="black", lty="dotted")
    rpy.r.legend("bottomleft", legend=[ str(bin) for bin in bins ], fill=cols, title="Bin size", horiz=True)
    rpy.r.dev_off()

    rpy.r.png("%s_variance_vs_major_allele.png" % unit, height=400, width=800)
    for bin, col in zip(profiles, cols)[0:-1]:
        profile = profiles[bin]
        plotfn = rpy.r.plot if bin == 0 else rpy.r.lines
        plotfn(y=[ na(tup[2]) for tup in profile ], x=[ tup[0] for tup in profile ], type='l', col=col,
            xlab="Major allele length (units)",
            ylab="Variance(distn of stutter)",
            main="%s: dependence of variance on major allele length" % unit)

    # bin=len, aggregate across all major allele lengths
    profile = profiles[max(profiles)]
    rpy.r.lines(y=[ na(tup[2]) for tup in profile ], x=[ tup[0] for tup in profile ], type='l', col="black", lty="dotted")
    rpy.r.legend("topleft", horiz=True, legend=[ str(bin) for bin in bins ], fill=cols, title="Bin size")
    rpy.r.dev_off()


import rpy
import sys
import simplejson
from pprint import pprint
from collections import defaultdict

if __name__ == "__main__":
    with open('several_singlecells/all_but_112082.model_params.json', 'r') as f:
        m = simplejson.load(f)

    # Look at how mean and variance of MSI depends on major allele length
    for unit in m['counts']:
        # The first bin (=0) means no binning, the final bin (max) means combine
        # data from all major allele lengths, which loses the dependency.
        bins = [ 0, 1, 2, 3, 4, 5, max(int(allele) for allele in m['counts'][unit]) ]
        profiles = build_profiles_for_major_alleles(m['counts'][unit], bins)
        pprint(dict(profiles))
        plot_mean_and_var_against_major_allele(profiles, unit)
