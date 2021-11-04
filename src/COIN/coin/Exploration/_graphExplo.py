# Data structures:
# n samples and N nodes.
# We need those input data files.
#  - One describing the population structure. It's a vector of length n with
#    categorical values.
#  - Another with the nodes by sample pattern of presence absence. It's an
#    N by n binary matrix.
#  - Another describing the phenotypes. It's a vector of length n with binary
#    values.
#  - A final one describing the graph. it is a file containing the edges. Choice
#    of structure is still unknown. For now we assume that it is a matrix with
#    two columns, where each line is an edge represented by the two nodes (one
#    per column) that the edge connects.


import sys, copy, os
from coin.Exploration._graph_Westbrook import *
from coin.Exploration._helper_functions import *
from caldera.Statistics._Pop import *
from caldera.Statistics._Envelope import *
import numpy as np
from progressbar import ProgressBar
from scipy.stats import chi2

class structure:
    """Contains the main initial files for the analysis:
        Slots
        Pop: The vector of population
        neighbors: the list of edges
        pattern: the matrix of presence / absence of nodes by samples
        Pheno: the binary vector of categories
        lengths: the size in bp of the nodes
        """
    def __init__(self, Pop, neighbours, nodes, pheno, lengths):
        self.Pop = Pop
        self.neighbours = neighbours
        self.pattern = nodes
        self.Pheno = pheno
        self.lengths = lengths
    
    def ns(self):
        # We have k groups identified through k-means
        # We returns two vectors of length k with the number of sample in each
        # group, split by phenotype.
        clusters = np.unique(self.Pop)
        clusters.sort()
        n1s = np.zeros(len(clusters))
        n2s = np.zeros(len(clusters))
        
        for clus in clusters:
            n1s[clus] = sum(self.Pheno[self.Pop == clus])
            n2s[clus] = sum(self.Pop == clus) - n1s[clus]
        
        return n1s, n2s
    
class testables:
    """Contains the current list of testable subgraphs:
        Slots
        subgraphs: The testable subgraphs
        envelopes: The envelopes of those subgraphs
        next: the location of the next subgraph to add
        k: the current value of k
        kmax: the limit to how many subgraphs we find
        alpha: the value of alpha for FWER control
        TH: the current testability threshold based on alpha and k
        """
    def __init__(self, alpha, nNodes, kmax):
        self.subgraphs = np.empty(2, dtype = graph)
        self.pvals = np.empty(2, dtype = int)
        self.next = 0
        self.k = 1
        self.kmax = kmax
        self.alpha = alpha
        self.TH = chi2.isf(alpha, 1)
        self.itemtable = {k: [] for k in range(nNodes)}
    
    def stop(self):
        return self.k == self.kmax
    
    def increment(self):
        # We use k +1 instead of k since we first add the subgraph to the list
        # and then test whether we have too many testable subgraphs
        return self.next == self.k + 1
    
    def append(self, S):
        self.subgraphs[self.next] = S
        self.pvals[self.next] = S.pval
        self.next = self.next + 1
    
    def clean(self, prune = True):
        self.k = self.k + 1
        if self.stop():
            return
        self.TH = chi2.isf(self.alpha / self.k, 1)
        if prune:
            keep = self.pvals >= self.TH
        else:
            keep = np.ones(self.pvals.shape, dtype=bool)
        testable_subgraphs = self.subgraphs[keep]
        testable_pvals = self.pvals[keep]
        self.subgraphs = np.empty(self.k + 1, dtype = graph)
        self.pvals = np.empty(self.k + 1, dtype = np.float64)
        self.subgraphs[0:sum(keep)] = testable_subgraphs
        self.pvals[0:sum(keep)] = testable_pvals
        self.next = sum(keep)
    
    def enum(self, S, G, n1s, n2s, Lmax, verbose = False, prune = True):
        """Enumerate all subgraphs containing S:

        Add the subgraph to R if it is testable. Stop if we reach the threshold.
        Otherwise, continue exploring the graph G. 

        Args:
            self: the list of testable subgraphs we want to extend
            S: the subgraph
            G: the graph to which the subgraph belongs.
            n1s: An vector of group sizes among the first phenotype.
            n2s: An vector of group sizes among the second phenotype.
            Lmax: maximum size of the subgraph
            verbose: a boolean, default to False.
            prune: Whether to prune. Used for internal testing

        Returns:
            Nothing (modify R in place)

        """
        # We first compute the minimal p-value.
        S.compute_envelope(G.Pop, n1s, n2s)
        # If we are below the threshold, we add the subgraph to our list.
        if  S.pval >=  self.TH or not prune:
            # we check if it is closed
            closure = (G.pattern[[int(ne) for ne in S.neighbours]] | S.ys) != S.ys
            closure = closure.sum(axis = 1) == 0
            closure = list(closure)
            closed = not any(closure)
            if closed:
                self.append(S)
        # We then check if we hit the stopping criterion for the envelope.
        if S.Env < self.TH and prune:
            return
        if self.increment():
            # We increase k by one, we clean R and we resume the search
            self.clean(prune = prune)
            if self.stop():
                return
            if verbose:
                if self.k % 1000 == 0:
                    print("k = " + str(self.k), flush = True)
        if S.length > Lmax:
            return
        # If there are no stop, we continue exploring the graph in depth.
        for j in S.neighbours:
            if self.stop():
                break
            S2 = copy.deepcopy(S)
            j = int(j)
            S2.add_node(j, G.pattern[j,:], neighbours = G.neighbours[j],
                        length = G.lengths[j])
            explored = check_itemtable(j, S2.ys, self.itemtable)
            if explored:
                next
            if S2.candidates(j):
#                print(list(S2.squareNodes.keys()))
                self.enum(S = S2, G = G, n1s = n1s, n2s = n2s, Lmax = Lmax,
                          verbose = verbose, prune = prune)
            else:
                next
        return

def explore(G, alpha, n1s, n2s, kmax, Lmax = 10 ** 7, verbose = False):
    """Compute the set R by incremental values of k:
    
    Args:
        G: the graph we want to explore.
        alpha: between 0 and 1. FWER we want to control.
        n1s: An vector of group sizes among the first phenotype.
        n2s: An vector of group sizes among the second phenotype.
        kmax: an integer, the maximum value of k we allow.
        Lmax: maximum size of the subgraph
        verbose: a boolean, default to false. Passed to enum
    
    Return: 
        The array R of testable subgraphs.
    """
    # Initiate the list of candidates.
    n = int(n1s.sum() + n2s.sum())
    nNodes = G.pattern.shape[0]
    R = testables(alpha, nNodes, kmax)
    pbar = ProgressBar()
    
    # Starting from all individual nodes, enumerate the subgraphs.
    # This will append R as we go on.
    # This could be parallelized, depending on how we pass R around.
    for s in pbar(np.arange(np.shape(G.pattern)[0])):
        if R.stop():
            break
        S = graph(n = n)
        explored = check_itemtable(s, G.pattern[s,:], R.itemtable)
        if explored:
            next
        S.add_node(s, G.pattern[s,:], neighbours = G.neighbours[s],
                   length = G.lengths[s])
        # Then we go on enumerating all pairs.
        # Note: a single node has no Vertex separator.
        R.enum(S, G, n1s, n2s, verbose = verbose, Lmax = Lmax)
        
    # We have finished and we return R
    return R.subgraphs, R.pvals, R.k, R.TH
