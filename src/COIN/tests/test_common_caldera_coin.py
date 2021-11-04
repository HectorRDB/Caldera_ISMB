import numpy as np
import caldera.data._toyData as toyDataset
from caldera.Exploration._ExploreDFS import explore as caldera_explore
from caldera.Exploration._Helper_Explo import structure as caldera_structure
from coin.Exploration._graphExplo import explore as coin_explore
from coin.Exploration._graphExplo import structure as coin_structure

def run_coin(Pop, neighbours, pattern, Pheno, Lengths):
    G = coin_structure(Pop, neighbours, pattern, Pheno, Lengths)
    n1s, n2s = G.ns()
    R, _, _, _ = coin_explore(G, 1, n1s, n2s, 10 ** 8)
    R = R[R != None]
    return(R)

def run_caldera(Pop, neighbours, pattern, Pheno, Lengths):
    G = caldera_structure(Pop, neighbours, pattern, Pheno, Lengths)
    R, _, _ = caldera_explore(G, 1, 10 ** 8)
    R = R[R != None]
    return(R)

def test_commmon_1():
    _, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset(n = 50, N = 200)
    Pop = np.zeros((pattern.shape[1],), dtype = np.int)
    Lengths = np.ones(pattern.shape[0])
    R_caldera = run_caldera(Pop, neighbours, pattern, Pheno, Lengths)
    R_coin = run_coin(Pop, neighbours, pattern, Pheno, Lengths)
    assert(R_caldera.shape == R_coin.shape), "COIN and CALDERA do not return the same objects"
    caldera_nodes=[sorted(list(S.nodes)) for S in R_caldera]
    coin_nodes=[sorted(list(S.squareNodes.keys()))for S in R_coin]
    for nodes in caldera_nodes:
        assert(nodes in coin_nodes), "Not all subgraphs are the same"

def test_commmon_2():
    _, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset(n = 30, N = 100)
    Pop = np.zeros((pattern.shape[1],), dtype = np.int)
    Lengths = np.ones(pattern.shape[0])
    R_caldera = run_caldera(Pop, neighbours, pattern, Pheno, Lengths)
    R_coin = run_coin(Pop, neighbours, pattern, Pheno, Lengths)
    assert(R_caldera.shape == R_coin.shape), "COIN and CALDERA do not return the same objects"
    caldera_nodes=[sorted(list(S.nodes)) for S in R_caldera]
    coin_nodes=[sorted(list(S.squareNodes.keys()))for S in R_coin]
    for nodes in caldera_nodes:
        assert(nodes in coin_nodes), "Not all subgraphs are the same"