import copy
import numpy as np
from caldera.Statistics._Envelope import Th

class light_graph:
    """A light graph is a subgraph of the De Bruijn graph with minimal info
        stored
    Slots:
        nodes: the list of current nodes.
        minp: the subgraph's minimal p-value
        pval: the subgraph's p-value
        ys: the subgraph pattern
    """
    
    def __init__(self, S):
        """Copying a light subgraph from a non-light subgraph
        
        Args:
            S: a non-light subgraph
        Return:
            An light version of the subgraph.
        """
        # Create the subgraph with one node
        self.nodes = copy.copy(S.nodes)
        self.minp = copy.copy(S.minp)
        self.ys = copy.copy(S.ys)
        self.pval = 0
        return
    
    def TH(self, Pop, Pheno):
        """ Compute the p-value of the subgraph (up to a chi-square transformation)"""
        clusters = np.unique(Pop)
        clusters.sort()
        n1s = np.zeros(len(clusters))
        n2s = np.zeros(len(clusters))
        a_s = np.zeros(len(clusters))
        xs = np.zeros(len(clusters))
        
        for clus in clusters:
            a_s[clus] = self.ys[(Pop==clus) & (Pheno==1)].sum()
            xs[clus] = self.ys[Pop == clus].sum()
            n1s[clus] = sum(Pheno[Pop == clus])
            n2s[clus] = sum(Pop == clus) - n1s[clus]
        
        return Th(a_s, xs, n1s, n2s)
