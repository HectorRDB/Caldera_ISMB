################################################################################
################################################################################

import numpy as np
import random

################################################################################
# Group of functions linked to closure enumeration
################################################################################
def generate_phenotype(n, prop = .5):
    phenotype = np.random.choice([0, 1], size=(n,), p=[prop, 1 - prop])
    return phenotype

def generate_pattern(N, phenotype, power = .95, density = 10):
    n = phenotype.shape[0]
    sig = int(N / density)
    null_pattern = np.random.choice([0, 1], size=(N - sig, n), p=[.5, .5])
    sig_pattern = [np.empty((n,)) for i in range(sig)]
    sig_pattern = np.array(sig_pattern, dtype = np.int)
    groups = np.split(np.array(range(sig)), int(sig / 10))
    for nodes in groups:
        same = np.random.choice([0, 1], size=(n,), p=[.05, .95])
        global_sig_pat = phenotype * same + (1 - phenotype) * (1-same)
        nodes_sig_pat = [global_sig_pat for i in range(nodes.shape[0])]
        nodes_sig_pat = np.array(nodes_sig_pat)
        split = np.random.choice([0, 1], size=nodes_sig_pat.shape,
                                 p=[.1, .90])
        nodes_sig_pat = nodes_sig_pat * split
        nodes_sig_pat[nodes_sig_pat.shape[0] - 1, ] = (
            ((1 - (nodes_sig_pat.sum(axis = 0) > 1)) *
            nodes_sig_pat[nodes_sig_pat.shape[0] - 1, ]) +
            nodes_sig_pat[nodes_sig_pat.shape[0] - 1, ])
        sig_pattern[list(nodes),] = nodes_sig_pat
    sig_pattern = np.array(sig_pattern, dtype = np.int)
    pattern = np.concatenate([null_pattern, sig_pattern])
    return pattern


def generate_random_graph(N, density = 10):
    edges = []
    sig = int(N / density)
    for i in range(N - sig):
        nei = random.sample(range(N - sig), 3)
        for n2 in nei:
            if not i == n2:
                edges.append([i, n2])
    groups = np.split(np.array(range(N - sig, N)), int(sig / 10))
    for nodes in groups:
        for i, node in enumerate(nodes):
            nei = random.sample(range(N - sig), 2)
            for n2 in nei:
                if not node == n2:
                    edges.append([i, n2])
            if not i == nodes.shape[0] - 1:
                edges.append([node, nodes[i + 1]])
            else:
                edges.append([node, nodes[0]])
    edges = np.array(edges)
    neighbours = dict()
    for i in range(N):
        neighbours[i] = []
    for i in range(edges.shape[0]):
        edge = edges[i, :]
        neighbours[edge[0]].append(edge[1])
        neighbours[edge[1]].append(edge[0])
    return edges, neighbours, groups

def generate_pop(n, C = 3):
    pop = np.random.choice(range(C), size=(n,), p= np.ones((C,)) / C)
    return pop
