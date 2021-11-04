################################################################################
################################################################################

import sys, os, argparse, copy, time
import numpy as np
sys.path.append(os.path.abspath("/Analysis/Simulations"))
from Generate_data import *
sys.path.append(os.path.abspath("src/FACS"))
import ExploreBFS as bfs


if __name__ == '__main__':
    random.seed(12)
    np.random.seed(32)
    output='/Output/Simulations/Imbalance/'
    alpha = 10 ** (-20)
    Lmax = 10 ** 5
    kmax = 10 ** 8
    n = 100
    N = 200
    times = np.zeros((7, 2))
    for i, ratio in enumerate([0.01, 0.02, 0.03, 0.04, 0.05, 0.1, .2]):
        Pheno = generate_phenotype(n, .5)
        pattern = generate_pattern(N, Pheno)
        edges, neighbours, groups = generate_random_graph(N, 10)
        Pop = np.zeros(Pheno.shape, dtype = np.int)
        Lengths = np.ones((N,))
        G = bfs.structure(Pop, neighbours, pattern, Pheno, Lengths)
        prop = 1 / (1 + ratio)
        Pop = np.zeros(Pheno.shape, dtype = np.int)
        ones = np.where(Pheno == 1)[0]
        Pop[ones[0:int(len(ones) * prop)]] = 1
        zeros = np.where(Pheno == 0)[0]
        Pop[zeros[0:int(len(zeros) * prop)]] = 1
        G.Pop = Pop
        f = open(output + 'partial_' + str(ratio) + '.txt','w')
        sys.stdout = f
        tick = time.time()
        R, Env, k0, TH = bfs.explore(G = G, alpha = alpha, kmax = kmax,
                                     Lmax = Lmax, verbose = True,
                                     threads = 1, old_Env = True)
        tock = time.time()
        times[i, 0] = tock - tick
        f.close()
        f = open(output + 'full_' + str(ratio) + '.txt','w')
        sys.stdout = f
        tick = time.time()
        R, Env, k0, TH = bfs.explore(G = G, alpha = alpha, kmax = kmax,
                                     Lmax = Lmax, verbose = True,
                                     threads = 1, old_Env = False)
        tock = time.time()
        times[i, 1] = tock - tick
        f.close()
    np.savetxt(output + 'times.txt', times, fmt = '%d')

                
                
                
