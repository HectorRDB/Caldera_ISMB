import numpy as np
import random
import caldera.Statistics._Envelope as Envelope
from statsmodels.stats.contingency_tables import StratifiedTable

# Setup
n = 15
K = 2
pheno = np.random.choice([0, 1], size=(n), p = [.7, .3])
pattern = np.random.choice([0, 1], size=(n), p = [.4, .6])
pop = np.random.choice(np.arange(K), size=(n), p = np.ones(K) / K)
# Normal way
tables = np.zeros((2, 2, K))
for k in range(K):
    phenoPop = pheno[np.where(pop == k)]
    patternPop = pattern[np.where(pop == k)]
    tables[0, 0, k] = len(np.where(phenoPop[np.where(patternPop == 1)] == 1)[0])
    tables[0, 1, k] = len(np.where(phenoPop[np.where(patternPop == 0)] == 1)[0])
    tables[1, 0, k] = len(np.where(phenoPop[np.where(patternPop == 1)] == 0)[0])
    tables[1, 1, k] = len(np.where(phenoPop[np.where(patternPop == 0)] == 0)[0])
# Our implementation
a = tables[0, 0, :]
xs = np.sum(tables, axis = 0)[0, :]
n1s = np.sum(tables, axis = 1)[0, :]
n2s = np.sum(tables, axis = 1)[1, :]

def test_cmd_test():
    pvalue = StratifiedTable(tables).test_null_odds().pvalue
    pvalue2 = Envelope.chi2.sf(Envelope.Th(a.sum(), xs, n1s, n2s), 1)
    assert((pvalue - pvalue2) / pvalue2 < 10**(-10))

def test_minp():
    minP = Envelope.chi2.sf(Envelope.minimal_p_value(xs, n1s, n2s), 1)
    a0s = range(int(max(0, xs[0] - n2s[0])),
            int(min(xs[0], n1s[0])) + 1)
    a1s = range(int(max(0, xs[1] - n2s[1])),
                int(min(xs[1], n1s[1])) + 1)
    minPbrute = 1
    tables2 = np.zeros((2, 2, K))
    for a0 in a0s:
        for a1 in a1s:
            tables2[0, 0, 0] = a0
            tables2[1, 0, 0] = xs[0] - a0
            tables2[0, 1, 0] = n1s[0] - a0
            tables2[1, 1, 0] = n2s[0] + a0 - xs[0]
            tables2[0, 0, 1] = a1
            tables2[1, 0, 1] = xs[1] - a1
            tables2[0, 1, 1] = n1s[1] - a1
            tables2[1, 1, 1] = n2s[1] + a1 - xs[1]
            pvalue = StratifiedTable(tables2).test_null_odds().pvalue
            if pvalue < minPbrute:
                minPbrute = pvalue
    assert(((minP - minPbrute) / minPbrute) < 10**(-10))

def test_envelope():
    n1s = np.array([4, 4])
    n2s = np.array([5, 5])
    Xs = np.array([6, 3])
    minP = Envelope.chi2.sf(Envelope.envelope(Xs, n1s, n2s), 1)
    ns = n1s + n2s
    x0s = range(Xs[0], ns[0] + 1)
    x1s = range(Xs[1], ns[1] + 1)
    minPbrute = 1
    for x0 in x0s:
        for x1 in x1s:
            xs = np.array([x0, x1])
            a0s = range(int(max(0, xs[0] - n2s[0])),
                        int(min(xs[0], n1s[0])) + 1)
            a1s = range(int(max(0, xs[1] - n2s[1])),
                        int(min(xs[1], n1s[1])) + 1) 
            tables2 = np.zeros((2, 2, K))
            for a0 in a0s:
                for a1 in a1s:
                    tables2[0, 0, 0] = a0
                    tables2[1, 0, 0] = xs[0] - a0
                    tables2[0, 1, 0] = n1s[0] - a0
                    tables2[1, 1, 0] = n2s[0] + a0 - xs[0]
                    tables2[0, 0, 1] = a1
                    tables2[1, 0, 1] = xs[1] - a1
                    tables2[0, 1, 1] = n1s[1] - a1
                    tables2[1, 1, 1] = n2s[1] + a1 - xs[1]
                    if not (tables2[0:2, 1, 0:2] == 0).all():
                        pvalue = StratifiedTable(tables2).test_null_odds().pvalue
                    if pvalue < minPbrute:
                        minPbrute = pvalue
    assert(((minP - minPbrute) / minPbrute) < 10**(-10)) 

