_, neighbours, pattern, Pheno, edges = toyDataset.generateToyDataset(n = 10, N = 50)
Pop = np.zeros((pattern.shape[1],), dtype = np.int)
Lengths = np.ones(pattern.shape[0])
G = structure(Pop, neighbours, pattern, Pheno, Lengths)

batch_size = 2
threads = 1
backend = "multiprocessing"
sMax = 5
verbose = False
Lmax = 10 ** 7
alpha=1
kmax=10**8
n1s, n2s = G.ns()
n = int(n1s.sum() + n2s.sum())
nNodes = G.lengths.shape[0]
temp_folder = G.offload(verbose = verbose)
from caldera.Exploration._ExploreBFS import *

sols = solutions(alpha, nNodes, Lmax, kmax, sMax)
with parallel_backend(backend, n_jobs=threads):
    C, stage = sols.build_first_stage(threads, G, n1s, n2s, verbose)

self = sols 
     
# [9, 11, 33, 47]