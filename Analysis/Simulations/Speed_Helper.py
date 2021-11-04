################################################################################
################################################################################

import sys, os, time, random, tracemalloc
import numpy as np
sys.path.append(os.path.abspath("/Analysis/Simulations/"))
from Generate_data import *
import caldera.Exploration._ExploreDFS as dfs
import caldera.Exploration._ExploreBFS as bfs
import caldera.Exploration._Helper_Explo as Helper_Explo
import coin.Exploration._graphExplo as graphExplo
from multiprocessing import Process

def BFS_Run(G, alpha, kmax, Lmax, threads, backend = "multiprocessing"):
    R, k0, TH = bfs.explore(G = G, alpha = alpha, kmax = kmax,
                                 Lmax = Lmax, verbose = True,
                                 threads = threads, backend = backend)
    return k0


def DFS_Run(G, alpha, kmax, Lmax):
    R, k0, TH = dfs.explore(G = G, alpha = alpha, kmax = kmax,
                             Lmax = Lmax, verbose = True)
    return k0


def COIN_Run(G, alpha, kmax, Lmax):
    n1s, n2s = G.ns()
    R, Env, k0, TH = graphExplo.explore(G = G, n1s = n1s, n2s = n2s,
                                        alpha = alpha, kmax = kmax,
                                        Lmax = Lmax, verbose = True)
    return k0
                             

def run_simu(id, n = 100, alpha = 0.05, p = 0.5, C = 1,
             timeout = 24 * 60 * 60, cutoff = 100):
    output='/Output/Simulations/Speed/'
    Lmax = 50 ** 5
    kmax = 10 ** 8
    res = np.zeros((4))
    for i, N in  enumerate([100, 200, 500, 1000, 2000, 5000, 10000, 20000]):
        print(N, flush = True)
        Pheno = generate_phenotype(n, p)
        pattern = generate_pattern(N, Pheno)
        edges, neighbours, groups = generate_random_graph(N, 10)
        Pop = generate_pop(n = n, C = C)
        Lengths = np.ones((N,))
        G = Helper_Explo.structure(Pop, neighbours, pattern, Pheno, Lengths)
        if N >= cutoff:
            start_time = time.time()
            BFS_Run(G, alpha, kmax, Lmax, threads = 1)
            duration = round(time.time() - start_time, 2)
            res[0] = duration
            print("Finished BFS. It took {t} seconds".format(t = duration))
            P = Process(target=DFS_Run, name='DFS', args = (G, alpha, kmax, Lmax))
            P.start()
            start_time = time.time()
            P.join(timeout = timeout)
            P.terminate()
            if P.exitcode == None:
                duration = timeout + 1
            elif P.exitcode < 0:
                duration = timeout + 1
            else:
                duration = time.time() - start_time
            res[1] = duration
            print("Finished DFS. It took {t} seconds".format(t = duration))
            P = Process(target=COIN_Run, name='COIN', args = (G, alpha, kmax, Lmax))
            P.start()
            start_time = time.time()
            P.join(timeout = timeout)
            P.terminate()
            if P.exitcode == None:
                duration = timeout + 1
            elif P.exitcode < 0:
                duration = timeout + 1
            else:
                duration = time.time() - start_time
            res[2] = duration
            print("Finished COIN. It took {t} seconds".format(t = duration))
            start_time = time.time()
            BFS_Run(G, alpha, kmax, Lmax, threads = 5)
            duration = round(time.time() - start_time, 2)
            res[3] = duration
            print("Finished BFS 5 threads. It took {t} seconds".format(t = duration))
            np.savetxt(output + "{N}_".format(N = N) + "res_" + str(id) + ".txt",
                       res, fmt = '%d')


def run_long_simu(id, n = 100, alpha = 0.05, p = 0.5, C = 1,
                  timeout = 24 * 60 * 60, cutoff = 100):
    output='/Output/Simulations/Speed/'
    Lmax = 50 ** 5
    kmax = 10 ** 8
    res = np.zeros((4))    
    N = 100000
    print(N, flush = True)
    Pheno = generate_phenotype(n, p)
    pattern = generate_pattern(N, Pheno)
    edges, neighbours, groups = generate_random_graph(N, 10)
    Pop = generate_pop(n = n, C = C)
    Lengths = np.ones((N,))
    G = Helper_Explo.structure(Pop, neighbours, pattern, Pheno, Lengths)
    start_time = time.time()
    BFS_Run(G, alpha, kmax, Lmax, threads = 5)
    duration = round(time.time() - start_time, 2)
    print(duration, flush = True)