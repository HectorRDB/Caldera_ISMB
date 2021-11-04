################################################################################

# The functions below are all used to efficiently explore the graph while 
# computing the fewest possible number of envelopes.

################################################################################
from coin.Exploration._graph_Westbrook import *
import caldera.Statistics._Envelope as Envelope
import numpy as np
import copy
from scipy.stats import chi2

def check_itemtable(n, pattern, itemtable):
    itemsets = itemtable[n]
    explored=False
    for explored_pat in itemsets:
        explored = all((explored_pat | pattern) == pattern)
        if explored:
            return explored
    itemtable[n].append(pattern)
    return explored

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
        min_p = Envelope.minimal_p_value(freqS, n1s, n2s)
    else:
        min_p = np.array([compute_min_ps([S], Pop, n1s, n2s) for S in Ss],
                         dtype = np.float64)
    return min_p
