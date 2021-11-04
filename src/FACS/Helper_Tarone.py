################################################################################
# The functions below are all used to efficiently explore the graph while 
# computing the fewest possible number of envelopes. 
################################################################################

import numpy as np
from scipy.stats import chi2
from Subgraphs import *

################################################################################
# Group of functions linked to Tarone's pruning
################################################################################
def compute_min_ps(Ss, Pop, n1s, n2s):
    """Compute the minimal p-value of a subgraph

    Args:
        Ss: an array of objects of class graph for which we want to compute the
        envelope
        Pop: the vector of populations
        n1s: An vector of group sizes among the first phenotype.
        n2s: An vector of group sizes among the second phenotype.

    Return:
        The minimal p-value of the each subgraph
    """
    if len(Ss) == 1:
        S = Ss[0]
        freqS = S.frequencies(Pop)
        min_p = minimal_p_value(freqS, n1s, n2s)
    else:
        min_p = np.array([compute_min_ps([S], Pop, n1s, n2s) for S in Ss],
                         dtype = np.float64)
    return min_p

def find_ko(min_ps, alpha, k = 1):
    """Compute the k0 according to the method in Terada et al
    
    Args:
        min_ps: The list of minimal p-values
        alpha: The threshold for controlling the FWER
        k: the previous value of k0 computed on a smaller set. Default to 1

    Return:
        The k0 value and associated cutoff for the Chi-square test
    """
    sortedMinP = -np.sort(-min_ps)
    N = min_ps.shape[0]
    a = k - 1
    b = N
    while b - a > 1:
        mid = int((a + b) / 2)
        TH = chi2.isf(alpha / (mid + 1), 1)
        if sortedMinP[mid] >= TH:
            a = mid 
        else:
            b = mid
    k0 = a + 1
    TH = chi2.isf(alpha / k0, 1)
    return k0, TH
    
