################################################################################
################################################################################

import sys, os, time, random, tracemalloc
import numpy as np
sys.path.append(os.path.abspath("/Analysis/Simulations/"))
from Generate_data import *
from Speed_Helper import BFS_Run, DFS_Run, COIN_Run
import caldera.Exploration._ExploreDFS as dfs
import caldera.Exploration._ExploreBFS as bfs
import caldera.Exploration._Helper_Explo as Helper_Explo
import coin.Exploration._graphExplo as graphExplo
from multiprocessing import Process
                             

if __name__ == '__main__':
    random.seed(12345)
    np.random.seed(32)
    output='/Output/Simulations/Speed/'
    alpha = 0.05
    Lmax = 50 ** 5
    kmax = 10 ** 8
    n = 50
    res = np.zeros((3))
    for i, N in  enumerate([100, 200, 500, 1000, 2000, 5000, 10000, 20000]):
        print(N, flush = True)
        Pheno = generate_phenotype(n, .5)
        pattern = generate_pattern(N, Pheno)
        edges, neighbours, groups = generate_random_graph(N, 10)
        Pop = np.zeros(Pheno.shape, dtype = np.int)
        Lengths = np.ones((N,))
        G = Helper_Explo.structure(Pop, neighbours, pattern, Pheno, Lengths)
        tracemalloc.start()
        BFS_Run(G, alpha, kmax, Lmax, threads = 1)
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        res[0] = peak / 10 ** 6
        print("Finished BFS. Peak memory usage was {m}G".format(m = peak / 10 ** 6))
        tracemalloc.start()
        DFS_Run(G, alpha, kmax, Lmax)
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        res[1] = peak / 10 ** 6
        print("Finished DFS. Peak memory usage was {m}G".format(m = peak / 10 ** 6))
        tracemalloc.start()
        COIN_Run(G, alpha, kmax, Lmax)
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        res[2] = peak / 10 ** 6
        print("Finished COIN. Peak memory usage was {m}G".format(m = peak / 10 ** 6))
        np.savetxt(output + "{N}_".format(N = N) + "mem.txt", res, fmt = '%d')
