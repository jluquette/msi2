from scipy.optimize import minimize_scalar
from utils import *

def get_weighted_distn(readhist, distn1, distn2):
    """Given probability distributions distn1 and distn2, return a function
    that combines the two distributions using a weight factor w."""
    # log-likelihood function--since many of the p-values will be very
    # close to 0, the log-likelihood will keep the sum much larger.  Very
    # small values could prove difficult for an optimizing algorithm to
    # handle.
    return lambda w: sum(nreads * log(w*distn1(allele) + (1-w)*distn2(allele))
                         for allele, nreads in readhist.iteritems())
    # The likelihood function
    #return lambda w: prod((w*f1(allele) + (1-w)*f2(allele))**nreads
                          #for allele, nreads in readhist.iteritems())


def genotype(readhist, pmfs, chisq_cutoff=0.99):
    """Genotype a locus with data in readhist given a set of probability mass
    functions `pmfs` which depends on major allele length (NOT repeat unit)."""

    pairs = combinations_with_replacement(readhist.keys()):
    hom_models = {}
    raw_results = []
    for allele1, allele2 in pairs:
        dist = get_weighted_distn(readhist, pmfs[allele1], pmfs[allele2])
        # Find the weight that produces the max likelihood for this genotype
        # (minimize the negative value)
        min_res = minimize_scalar(lambda x: -dist(x), bounds=(0,1), method='bounded')
        # (likelihood, weight, allele1 len, allele2 len)
        res = (-min_res.fun, min_res.x, allele1, allele2)
        raw_results.append(res)
        if allele1 == allele2:
            hom_models[allele1] = res
    print(raw_results)

    # If allele1 = allele2 (a homozygous model), then the weight parameter
    # in the weighted distn drops out.  This means homozygous models are at
    # a disadvantage compared against heterozygous models since the het models
    # have a free weight parameter that is estimated from the read data.  We
    # apply a log-ratio test to see if the added parameter is justified.
    adjusted_results = []
    for likelihood, weight, allele1, allele2 in raw_results:
        # Only keep the 2 allele model if it's significantly better than both
        # of the single allele models.
        l1 = hom_models[allele1][0]
        stat1 = -2 * log(l1/likelihood)
        p1 = rpy.r.pchisq(stat1, df=1)

        l2 = hom_models[allele2]
        stat2 = -2 * log(l2/likelihood)
        p2 = rpy.r.pchisq(stat2, df=1)

        # The larger the chi square statistic, the smaller the single allele's
        # likelihood was compared to the double allele model.
        # We don't get away with cases where allele1=allele2.  In those cases,
        # l1 = l2 = likelihood, so we're taking log(1)=0, so our statistic
        # is 0.  But the chi square probability mass function 
        if p1 < chisq_cutoff and p2 < chisq_cutoff or allele1 == allele2:
            adjusted_results.append((likelihood, weight, allele1, allele2))

    # Get the best one
    return sorted(adjusted_results)[0]


if __name__ == "__main__":
    import simplejson

    with open(sys.args[1], 'r') as f:
        data = simplejson.load(f)['counts']

    mono_data = {}
    for allele in data['mono']:
        mono_data[allele] = data['mono'][allele]['Hist']

    mono_pmfs = {}
    mono_major_alleles = [ int(k) for k in mono_pmfs.keys() ]
    for allele in range(min(mono_major_alleles), max(mono_major_alleles)+ 1):
        mono_pmfs[allele] = 

    test_data = [
        (13, { '0': 15 }),          # Obvious homozygous
        (13, { '0': 15, '1': 2 }),  # Make sure the chisq test rejects het here
        (13, { '-2': 1, '-1': 5, '3': 10, '4': 5 })
    ]
