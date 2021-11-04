import numpy as np
import random
def generateToyDataset(seed = 571, n = 100, N = 500):
    np.random.seed(seed)
    random.seed(a = seed)
    alpha = .05
    # Simulate the phenotype
    Pheno = [0] * n + [1] * n
    np.random.shuffle(Pheno)
    Pheno = np.array(Pheno)
    # Simulate the pop structure
    Pop = []
    for i in range(0, 4):
        Pop = Pop + [i] * 50
    np.random.shuffle(Pop)
    Pop = np.array(Pop)
    # Simulate the pattern
    pattern = np.random.choice([0, 1], size=(N,2 * n), p=[.5, .5])
    # Simulate the edges
    edges = []
    for i in range(2 * N):
        n1, n2 = random.sample(range(N), 2)
        if not n1 == n2:
            edges.append([n1, n2])
    edges = np.array(edges)
    edges = np.unique(edges, axis = 1)
    neighbours = dict()
    for i in range(N):
        neighbours[i] = []
    for i in range(2 * N):
        edge = edges[i, :]
        neighbours[edge[0]].append(edge[1])
        neighbours[edge[1]].append(edge[0])
    return Pop, neighbours, pattern, Pheno, edges