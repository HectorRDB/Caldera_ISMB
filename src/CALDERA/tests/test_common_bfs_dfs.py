import numpy as np
import caldera.data._toyData as toyDataset
from caldera.Exploration._ExploreDFS import explore as dfs_explore
from caldera.Exploration._ExploreBFS import explore as bfs_explore
from caldera.Exploration._Helper_Explo import *
from caldera.Statistics._Helper_Tarone import compute_TH

def run_dfs(G):
    R, k0, TH = dfs_explore(G, 1, 10 ** 8)
    R = R[R != None]
    return(R)

def run_bfs(G, batch_size = 10 ** 10):
    R, k0, TH = bfs_explore(G, 1, 10 ** 8, sMax = 50, batch_size = batch_size)
    R = R[R != None]
    return(R)

def test_commmon_1():
    _, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset(n = 10, N = 50)
    Pop = np.zeros((pattern.shape[1],), dtype = int)
    Lengths = np.ones(pattern.shape[0])
    G = structure(Pop, neighbours, pattern, Pheno, Lengths)
    R_DFS = run_dfs(G)
    R_BFS = run_bfs(G)
    R_BFS_Batch = run_bfs(G, batch_size=2)
    assert(R_DFS.shape == R_BFS.shape), "Both search returns the same objects"
    assert(R_DFS.shape == R_BFS_Batch.shape), "Both search returns the same objects"
    bfs_nodes = [sorted(list(S.nodes)) for S in R_BFS]
    bfs_batch_nodes = [sorted(list(S.nodes)) for S in R_BFS_Batch]
    dfs_nodes=[sorted(list(S.nodes)) for S in R_DFS]
    for nodes in bfs_nodes:
        assert(nodes in dfs_nodes), "Not all subgraphs are the same"
        assert(nodes in bfs_batch_nodes), "Not all subgraphs are the same"

def test_commmon_2():
    _, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset(n = 35, N = 200)
    Pop = np.zeros((pattern.shape[1],), dtype = int)
    Lengths = np.ones(pattern.shape[0])
    G = structure(Pop, neighbours, pattern, Pheno, Lengths)
    R_BFS_Batch = run_bfs(G, batch_size=1)
    R_DFS = run_dfs(G)
    R_BFS = run_bfs(G)
    assert(R_DFS.shape == R_BFS.shape), "Both search returns the same objects"
    assert(R_DFS.shape == R_BFS_Batch.shape), "Both search returns the same objects"
    bfs_nodes = [sorted(list(S.nodes)) for S in R_BFS]
    bfs_batch_nodes = [sorted(list(S.nodes)) for S in R_BFS_Batch]
    dfs_nodes=[sorted(list(S.nodes)) for S in R_DFS]
    for nodes in bfs_nodes:
        assert(nodes in dfs_nodes), "Not all subgraphs are the same"
        assert(nodes in bfs_batch_nodes), "Not all subgraphs are the same"


def test_common_p_value_comput():
    _, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset(n = 35, N = 200)
    Pop = np.zeros((pattern.shape[1],), dtype = int)
    Lengths = np.ones(pattern.shape[0])
    G = structure(Pop, neighbours, pattern, Pheno, Lengths)
    R = run_bfs(G)
    pvals = np.array([S.TH(G.Pop, G.Pheno) for S in R])
    # pvals = pvals[:,0]
    pvals2 = compute_TH(R, G.Pop, G.Pheno)
    pvals2 = [S.pval for S in pvals2]
    assert((pvals == pvals2).all()), 'the pvalues computation are wrong'


