################################################################################
# Implementation of our algorithm to explore all connected closed
# subgraphs in depth
################################################################################

import sys
import numpy as np
from scipy.stats import chi2
from caldera.Statistics._Envelope import *
from caldera.Statistics._Helper_Tarone import *
from caldera.Exploration._Helper_Explo import *
from caldera.Exploration._Subgraphs import *
from caldera.Exploration._Light_Subgraphs import *

################################################################################
# The solution stores all the info we need on the current set of testables
# subgraphs: the value of k, the subgraphs, their p-values and the fixed stoping
# points
################################################################################

class solutions:
    """Used to store the solutions as we build them:
        Slots
        alpha: The values of alpha
        k: the current value of k0
        TH: the corresponding chi-square
        R: the array of testable subgraphs
        Env: their envelopes
        Lmax: maximum allowed size of subgaph
        kmax: maximum allowed value of k
        """
        
    def __init__(self, alpha, nNodes, Lmax, kmax):
        self.alpha = alpha
        self.Lmax = Lmax
        self.kmax = kmax
        self.k = 1
        self.TH = chi2.isf(alpha, 1)
        self.R = np.array([], dtype = graph)
        self.Minps = np.array([], dtype = np.float64)
        
    def enum(self, S, G, verbose, n1s, n2s):
        # We prune if k is too big
        if self.k > self.kmax:
            return
        # We add the k-testable subgraph
        if S.minp >= self.TH:
            self.R = np.append(self.R, [light_graph(S)])
            self.Minps = np.append(self.Minps, [S.minp])
        # We prune if it is prunable and not-k-testable
        if S.Env < self.TH:
            return
        # If we need to, we update k0
        if len(self.R) > self.k:
            self.k = self.k + 1
            if verbose:
                if self.k % 1000 == 0:
                    print("k = " + str(self.k), flush = True)
                
            self.TH = chi2.isf(self.alpha / self.k, 1)
            testable = self.Minps >= self.TH
            self.R = self.R[testable]
            self.Minps = self.Minps[testable]
        # If the subgraph is too big, we stop
        if S.length > self.Lmax:
            return
        # We iterate for the parents
        for S2 in S.Childrens(G.Pop, G.pattern, G.neighbours, G.lengths,
                              self.TH, n1s, n2s):
            self.enum(S2, G, verbose, n1s, n2s)

################################################################################
# This is the part where we explore the graph per say
################################################################################
def explore(G, alpha, kmax, Lmax = 10 ** 7, verbose = False):
    """Compute the set R by incremental values of k:
    
    Args:
        G: the graph we want to explore.
        alpha: between 0 and 1. FWER we want to control.
        n1s: An vector of group sizes among the first phenotype.
        n2s: An vector of group sizes among the second phenotype.
        kmax: an integer, the maximum value of k we allow.
        Lmax: maximum size of the subgraph
        verbose: a boolean, default to false.
        threads: Number of cores used in the process
        output: where to store the output
    
    Return: 
        The array R of testable subgraphs.
    """
    
    # We avoid computing this too often
    n1s, n2s = G.ns()
    n = int(n1s.sum() + n2s.sum())
    nNodes = G.lengths.shape[0]
    
    # Initiate the list of candidates. We create the closure of every subgraph
    # of size 1
    sols = solutions(alpha, nNodes, Lmax, kmax)
    candidates = start(range(nNodes), G.Pop, G.pattern, G.neighbours, G.lengths,
                       sols.TH, n1s, n2s)
                       
    for S in candidates:
        sols.enum(S, G, verbose, n1s, n2s)
        
    # We have finished and we return the testable graphs
    return sols.R, sols.k, sols.TH
