from coin.Exploration._graphExplo import *
import caldera.data._Ground_truth as ground_truth

def run_dfs_explo(G, TH, n1s, n2s):
    R, _, _, _ = explore(G, 1, n1s, n2s, 10 ** 8)
    R = R[R != None]
    return(R)

def test_bfs_data_1():
    G, TH, n1s, n2s = ground_truth.generateToyData1(BFS = False)
    R = run_dfs_explo(G, TH, n1s, n2s)
    assert(R.shape[0] == 10), "DFS Exploration of the dataset 1 is not as expected"

def test_bfs_data_2():
    G, TH, n1s, n2s = ground_truth.generateToyData2(BFS = False)
    R = run_dfs_explo(G, TH, n1s, n2s)
    assert(R.shape[0] == 6), "DFS Exploration of the dataset 2 is not as expected"

def test_bfs_data_3():
    G, TH, n1s, n2s = ground_truth.generateToyData3(BFS = False)
    R = run_dfs_explo(G, TH, n1s, n2s)
    assert(R.shape[0] == 11), "DFS Exploration of the dataset 3 is not as expected"

def test_bfs_data_4():
    G, TH, n1s, n2s = ground_truth.generateToyData4(BFS = False)
    R = run_dfs_explo(G, TH, n1s, n2s)
    assert(R.shape[0] == 13), "DFS Exploration of the dataset 4 is not as expected"
