import numpy as np
import pandas as pd
import random
from caldera.Exploration._Subgraphs import *
from caldera.Exploration._ExploreBFS import *
import caldera.data._toyData as toyDataset

Pop, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset()
alpha = .05
Lengths = np.ones(pattern.shape[0])
G = structure(Pop, neighbours, pattern, Pheno, Lengths)
n1s, n2s = G.ns()
TH = chi2.isf(.05 / 1, 1)

def test_subgraph_construction():
    for node in range(G.pattern.shape[0]):
        S = graph()
        S.new_graph(node, G.Pop, G.pattern, G.neighbours, G.lengths, TH, n1s, n2s)
        assert(np.array_equal(S.ys, G.pattern[node]))
        if len(S.nodes) == 1:
            assert(np.array_equal(S.common, G.pattern[node]))
        if len(G.neighbours[node]) > 0:
            G.pattern[G.neighbours[node][0],] = G.pattern[node]
            if G.neighbours[node][0] > node:
                node = G.neighbours[node][0]
            S = create_Graph(node, G.Pop, G.pattern, G.neighbours, G.lengths, TH, n1s, n2s)
            assert(S.Env == envelope(S.frequencies(G.Pop), n1s, n2s))
            assert(np.array_equal(S.ys, G.pattern[node]))
            assert(len(S.nodes) > 1)

def test_subgraph_children():
    for node in range(G.pattern.shape[0]):
        if len(G.neighbours[node]) > 0:
            S = create_Graph(node, G.Pop, G.pattern, G.neighbours, G.lengths, TH, n1s, n2s)
            if S.closed:
                C = S.Childrens(G.Pop, G.pattern, G.neighbours, G.lengths, TH, n1s, n2s)
                assert (not any([S2.iS == S.iS for S2 in C]))
                assert(all([S2.Env == envelope(S2.frequencies(G.Pop), n1s, n2s) for S2 in C]))

def test_pruning():
    node = 123
    G.pattern[node] = np.ones(G.pattern[node].shape)
    S = graph()
    S.new_graph(node, G.Pop, G.pattern, G.neighbours, G.lengths, TH, n1s, n2s)
    assert(S.Env == 0)
    assert(S.length == 0)
    assert(S.Env == envelope(S.frequencies(G.Pop), n1s, n2s))
    node = 1
    G.pattern[G.neighbours[node]] = np.ones(G.pattern[G.neighbours[node]].shape)
    S = graph()
    S.new_graph(node, G.Pop, G.pattern, G.neighbours, G.lengths, TH, n1s, n2s)
    assert(S.Childrens(G.Pop, G.pattern, G.neighbours, G.lengths, TH, n1s, n2s).shape == (0, ))
    assert(S.Env == envelope(S.frequencies(G.Pop), n1s, n2s))
