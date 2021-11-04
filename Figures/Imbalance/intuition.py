import os, sys
import numpy as np
sys.path.append(os.path.abspath("CALDERA/scripts"))
from Envelope import *
from scipy.stats import chi2

n1 = 260
n2 = 20
ns = np.array([n1, n2], dtype = np.int)
n1s = np.array([int(n1 / 2), int(n2/2)])
n2s = np.array([int(n1 / 2), int(n2/2)])
max_val = envelope(n1s, n1s, n2s)

RES = np.ones((int((n1 + 1) * (n2 + 1)), 4), dtype = np.float64)

for i in range(n1 + 1):
    for j in range(n2 + 1):
        xs = np.array([i, j])
        res = [i, j, envelope(xs, n1s, n2s), max_val]
        RES[i + j * (n1 + 1),] = res
        if i >= (n1 / 2) and j >= (n2 / 2):
            RES[i + j * (n1 + 1), 3] = old_envelope(xs, n1s, n2s)

np.savetxt('Figures/Imbalance/env_1.txt', RES, fmt = '%d')

n1 = 140
n2 = 140
ns = np.array([n1, n2], dtype = np.int)
n1s = np.array([int(n1 / 2), int(n2/2)])
n2s = np.array([int(n1 / 2), int(n2/2)])
max_val = envelope(n1s, n1s, n2s)

RES = np.ones((int((n1 + 1) * (n2 + 1)), 4), dtype = np.float64)

for i in range(n1 + 1):
    for j in range(n2 + 1):
        xs = np.array([i, j])
        res = [i, j, envelope(xs, n1s, n2s), max_val]
        RES[i + j * (n1 + 1),] = res
        if i >= (n1 / 2) and j >= (n2 / 2):
            RES[i + j * (n1 + 1), 3] = old_envelope(xs, n1s, n2s)

np.savetxt('Figures/Imbalance/env_2.txt', RES, fmt = '%d')
