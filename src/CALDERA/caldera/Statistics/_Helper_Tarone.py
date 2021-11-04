################################################################################
# The functions below are all used to efficiently explore the graph while 
# computing the fewest possible number of envelopes. 
################################################################################

import numpy as np
from scipy.stats import chi2
from caldera.Exploration._Subgraphs import *

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
    # In case of ties
    if a + 1 >= N:
        k0 = a + 1
    elif sortedMinP[a] == sortedMinP[b]:
        k0 = b + 1
    else:
        k0 = a + 1
    TH = chi2.isf(alpha / k0, 1)
    # 
    return k0, TH

def find_alpha(C, Pop, Pheno, top = 10):
    """Find alpha such that top subgraphs of C are significant
    
    Args:
        C: The list of subgraphs
        Pop: the vector of populations
        Pheno: The vector of phenotypes
        top: How many subgraphs should be significant

    Return:
        The value of alpha
    """
    pvals = -np.sort(-np.array([S.TH(Pop, Pheno) for S in C]))
    TH = pvals[top - 1]
    return min(1, chi2.sf(TH, 1) * len(C))
    

def compute_TH(Ss, Pop, Pheno, TH = 0):
    """Compute the p-value of a subgraph

    Args:
        Ss: an array of objects of class graph for which we want to compute the
        envelope
        Pop: the vector of populations
        Pheno: The vector of phenotypes
        TH: The cutoff value. Default to 0, i.e. no filtering
    """
    # Prepare the matrices
    clusters = np.unique(Pop)
    clusters.sort()
    n1s = np.zeros(len(clusters))
    n2s = np.zeros(len(clusters))
    for clus in clusters:
        n1s[clus] = sum(Pheno[Pop == clus])
        n2s[clus] = sum(Pop == clus) - n1s[clus]
    ns = n1s + n2s
    
    a_s = np.zeros((len(Ss), len(clusters)))
    xs = np.zeros((len(Ss), len(clusters)))
    Ys = np.array([S.ys for S in Ss])
    Ys_Pheno_p = np.multiply(Ys, Pheno)
    for clus in clusters:
        a_s_clus = np.multiply(Ys_Pheno_p, Pop == clus)
        a_s[:,clus] = np.sum(a_s_clus, axis = 1)
        xs_clus = np.multiply(Ys, Pop==clus)
        xs[:,clus] = np.sum(xs_clus, axis = 1)
    
    # Compute the p-values
    nss = (ns * (ns - 1))
    num = np.multiply(xs, n1s)
    num = np.divide(num, ns)
    num = np.sum(a_s, axis = 1) - np.sum(num, axis = 1)
    num = np.square(num)
    denum1 = np.multiply(n1s, n2s)
    denum1 = np.multiply(denum1, xs)
    denum2 = 1 - np.divide(xs, ns)
    denum = np.multiply(denum1, denum2)
    denum = np.divide(denum, nss)
    denum = np.sum(denum, axis = 1)
    pvals = np.divide(num, denum)

    # Filter and return the subgraphs
    Ss = Ss[pvals >= TH]
    pvals = pvals[pvals >= TH]
    for (pval, S) in zip(pvals, Ss):  
        S.pval = pval
    
    return Ss

