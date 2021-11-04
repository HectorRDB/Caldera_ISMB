# # Testing bfs_exploration
from caldera.Exploration._Subgraphs import *
from caldera.Exploration._ExploreBFS import *
import caldera.data._Ground_truth as ground_truth
import caldera.data._toyData as toyDataset

def run(G, batch_size = 10 ** 10, sMax = 50):
    R, _, _ = explore(G, 1, 10 ** 8, sMax = sMax, batch_size = batch_size)
    R = R[R != None]
    return R

def explo(Pop, Patterns, Neighbours, Lengths, n1s, n2s):
    nNodes = Patterns.shape[0]
    C = start(np.array(range(nNodes)), Pop, Patterns, Neighbours, Lengths, 0, n1s, n2s)
    R = C
    while C.shape[0] > 0:
        C = expand(C, Pop, Patterns, Neighbours, Lengths, 0, n1s, n2s)
        if C.shape[0] == 0:
            break
        R = np.concatenate([R, C])
    
    return R

# Testing that the epxloration in breadth works
def test_bfs_data_1():
    G, _, _, _ = ground_truth.generateToyData1()
    R = run(G)
    assert(R.shape[0] == 10), "BFS Exploration of the dataset 1 is not as expected"
    R = run(G, batch_size = 4)
    assert(R.shape[0] == 10), "BFS Exploration of the dataset 1 is not as expected"

def test_bfs_data_2():
    G, _, _, _ = ground_truth.generateToyData2()
    R = run(G)
    assert(R.shape[0] == 6), "BFS Exploration of the dataset 2 is not as expected"
    R = run(G, batch_size = 4)
    assert(R.shape[0] == 6), "BFS Exploration of the dataset 2 is not as expected"

def test_bfs_data_3():
    G, _, _, _ = ground_truth.generateToyData3()
    R = run(G)
    assert(R.shape[0] == 11), "BFS Exploration of the dataset 3 is not as expected"
    R = run(G, batch_size = 4)
    assert(R.shape[0] == 11), "BFS Exploration of the dataset 3 is not as expected"

def test_bfs_data_4():
    G, _, _, _ = ground_truth.generateToyData4()
    R = run(G)
    assert(R.shape[0] == 13), "BFS Exploration of the dataset 4 is not as expected"
    R = run(G, batch_size = 4)
    
# Testing that all options work
def test_closure_1():
    _, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset(n = 2, N = 10)
    Pop = np.zeros((pattern.shape[1],), dtype = int)
    G, _, n1s, n2s = ground_truth.format(Pop, neighbours, pattern, Pheno, edges)
    R = explo(G.Pop, G.pattern, G.neighbours, G.lengths, n1s, n2s)
    for S in R:
        closure = (pattern[S.neighbours,] | S.ys) == S.ys
        closure = closure.all(axis = -1)
        assert(not closure.any()), "BFS exploration returns non-closed graphs"

def test_closure_2():
    _, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset(n = 50, N = 200)
    Pop = np.zeros((pattern.shape[1],), dtype = int)
    G, _, n1s, n2s = ground_truth.format(Pop, neighbours, pattern, Pheno, edges)
    R = explo(G.Pop, G.pattern, G.neighbours, G.lengths, n1s, n2s)
    for S in R:
        closure = (pattern[S.neighbours,] | S.ys) == S.ys
        closure = closure.all(axis = -1)
        assert(not closure.any()), "BFS exploration returns non-closed graphs"

def test_depth_limit():
    _, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset(n = 20, N = 100)
    Pop = np.zeros((pattern.shape[1],), dtype = int)
    G, _, _, _ = ground_truth.format(Pop, neighbours, pattern, Pheno, edges)
    R = run(G)
    ref_nodes = [sorted(list(S.nodes)) for S in R]
    for s in range(1, 21):
        R_bfs = run(G, sMax = s)
        R_batch = run(G, batch_size = 2, sMax = s)
        assert(R.shape[0] >= R_bfs.shape[0]), "Exploration does stop early"
        assert(R_bfs.shape == R_batch.shape), "Both search returns the same objects"
        bfs_batch_nodes = [sorted(list(S.nodes)) for S in R_batch]
        bfs_nodes=[sorted(list(S.nodes)) for S in R_bfs]
        for nodes in bfs_nodes:
            assert(nodes in ref_nodes), "Not all subgraphs are the same"
            assert(nodes in bfs_batch_nodes), "Not all subgraphs are the same"

def test_no_alpha():
    G, _, _, _ = ground_truth.generateToyData1()
    R, _, _ = explore(G, None, 10 ** 8, sMax = 50)
    assert(R.shape[0] == 10), "BFS Exploration of the dataset 1 is not as expected"
    G, _, _, _ = ground_truth.generateToyData2()
    R, _, _ = explore(G, None, 10 ** 8, sMax = 50)
    assert(R.shape[0] == 6), "BFS Exploration of the dataset 2 is not as expected"
    G, _, _, _ = ground_truth.generateToyData3()
    R, _, _ = explore(G, None, 10 ** 8, sMax = 50)
    assert(R.shape[0] == 11), "BFS Exploration of the dataset 3 is not as expected"
    G, _, _, _ = ground_truth.generateToyData4()
    R, _, _ = explore(G, None, 10 ** 8, sMax = 50)
    assert(R.shape[0] == 13), "BFS Exploration of the dataset 4 is not as expected"

