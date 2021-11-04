from caldera.Exploration._Subgraphs import *
from caldera.Exploration._ExploreDFS import *
import caldera.data._Ground_truth as ground_truth
import caldera.data._toyData as toyDataset

def run_dfs_explo(G, TH, n1s, n2s):
    nNodes = G.lengths.shape[0]
    sols = solutions(1, nNodes, Lmax = 500, kmax = 100000)
    candidates = start(range(nNodes), G.Pop, G.pattern, G.neighbours, G.lengths,
                   sols.TH, n1s, n2s)
    for S in candidates:
        sols.enum(S, G, True, n1s, n2s)
    return(sols.R)

def test_dfs_data_1():
    G, TH, n1s, n2s = ground_truth.generateToyData1()
    R = run_dfs_explo(G, TH, n1s, n2s)
    assert(R.shape[0] == 10), "DFS Exploration of the dataset 1 is not as expected"

def test_dfs_data_2():
    G, TH, n1s, n2s = ground_truth.generateToyData2()
    R = run_dfs_explo(G, TH, n1s, n2s)
    assert(R.shape[0] == 6), "DFS Exploration of the dataset 2 is not as expected"

def test_dfs_data_3():
    G, TH, n1s, n2s = ground_truth.generateToyData3()
    R = run_dfs_explo(G, TH, n1s, n2s)
    assert(R.shape[0] == 11), "DFS Exploration of the dataset 3 is not as expected"

def test_dfs_data_4():
    G, TH, n1s, n2s = ground_truth.generateToyData4()
    R = run_dfs_explo(G, TH, n1s, n2s)
    assert(R.shape[0] == 13), "DFS Exploration of the dataset 4 is not as expected"
