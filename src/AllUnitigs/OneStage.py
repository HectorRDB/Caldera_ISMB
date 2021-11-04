import numpy as np
import pandas as pd
from caldera.Statistics._Envelope import *
from caldera.Statistics._Helper_Tarone import *
from caldera._core import create_G
from caldera.Exploration._Helper_Explo import create_Graph
from caldera.Exploration._Light_Subgraphs import light_graph
from caldera.Post._dbgwas_format import create_files


def test_all_unitigs(loc, output, communities, comFile, alpha, verbose):
    G = create_G(loc, communities, output, comFile, verbose, restart = False)
    n1s, n2s = G.ns()
    Ss = [create_Graph(s, G.Pop, G.pattern, G.neighbours, G.lengths, -1, n1s, n2s)
         for s in range(G.lengths.shape[0])]
    Minps = np.array([S.minp for S in Ss])
    k0, TH = find_ko(Minps, alpha)
    are_testable = Minps >= TH
    Ss = np.array(Ss)[are_testable]
    R = np.array([light_graph(S) for S in Ss])
    pvals = np.array([S.TH(G.Pop, G.Pheno) for S in R])
    sig_graphs = pvals >= TH
    if np.sum(sig_graphs) == 0:
        print("We found no significant subgraphs", flush = True)
        return
    else:
        sigs = R[sig_graphs]
        sigs_pvals = pvals[sig_graphs]
    create_files(loc, output, sigs, sigs_pvals, G, verbose)
    