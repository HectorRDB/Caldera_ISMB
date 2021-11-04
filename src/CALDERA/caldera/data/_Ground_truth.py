import numpy as np
import random
from caldera.Exploration._Helper_Explo import structure

def format(Pop, neighbours, pattern, Pheno, edges):
    Lengths = np.ones(pattern.shape[0])
    G = structure(Pop, neighbours, pattern, Pheno, Lengths)
    n1s, n2s = G.ns()
    return G, 0, n1s, n2s

def generateToyData1():
    Pheno = np.array([0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0],
                     dtype = int)
    # Simulate the phenotype
    Pop = np.zeros((12,), dtype = int)
    pattern = np.array([[1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]])

    # Simulate the edges
    edges = [[0, 1],[1, 2], [2, 3], [3, 4], [4, 5],
             [5, 6]]
    neighbours = dict()

    for i in range(7):
        neighbours[i] = []

    for i in range(6):
        edge = edges[i]
        neighbours[edge[0]] = neighbours[edge[0]] + [edge[1]]
        neighbours[edge[1]] = neighbours[edge[1]] + [edge[0]]

    return format(Pop, neighbours, pattern, Pheno, edges)


def generateToyData2():
    Pheno = np.array([0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0],
                     dtype = int)
    # Simulate the phenotype
    Pop = np.zeros((12,), dtype = int)
    pattern = np.array([[0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]])

    # Simulate the edges
    edges = [[0, 1],[1, 2]]
    neighbours = dict()

    for i in range(3):
        neighbours[i] = []

    for i in range(2):
        edge = edges[i]
        neighbours[edge[0]] = neighbours[edge[0]] + [edge[1]]
        neighbours[edge[1]] = neighbours[edge[1]] + [edge[0]]

    return format(Pop, neighbours, pattern, Pheno, edges)

def generateToyData3():
    Pheno = np.array([0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0],
                     dtype = int)
    # Simulate the phenotype
    Pop = np.zeros((12,), dtype = int)
    pattern = np.array([[0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1]])

    # Simulate the edges
    edges = [[0, 2],[1, 2], [2, 3]]
    neighbours = dict()

    for i in range(4):
        neighbours[i] = []

    for i in range(3):
        edge = edges[i]
        neighbours[edge[0]] = neighbours[edge[0]] + [edge[1]]
        neighbours[edge[1]] = neighbours[edge[1]] + [edge[0]]

    return format(Pop, neighbours, pattern, Pheno, edges)

def generateToyData4():
    Pheno = np.array([0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0],
                     dtype = int)
    # Simulate the phenotype
    Pop = np.zeros((12,), dtype = int)
    pattern = np.array([[0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                        [0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]])

    # Simulate the edges
    edges = [[0, 2],[1, 2], [2, 3], [1, 4]]
    neighbours = dict()

    for i in range(5):
        neighbours[i] = []

    for i in range(4):
        edge = edges[i]
        neighbours[edge[0]] = neighbours[edge[0]] + [edge[1]]
        neighbours[edge[1]] = neighbours[edge[1]] + [edge[0]]

    return format(Pop, neighbours, pattern, Pheno, edges)
