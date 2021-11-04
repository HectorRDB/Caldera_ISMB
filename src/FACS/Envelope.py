# This provide all the necessary functions to compute the minimal p-value and
# envelope

import numpy as np
from scipy.stats import chi2

def Th(a, xs, n1s, n2s):
    ns = n1s + n2s
    num = np.multiply(xs, n1s)
    num = np.divide(num, ns)
    num = a - num.sum()
    num = num**2
    denum = n1s * n2s * xs * (1 - xs / ns) / (ns * (ns - 1))
    denum = denum.sum()
    
    return num / denum
    
def minimal_p_value(xs, n1s, n2s):
    """Given the frequencies in the various populations, compute the minimal p-value."""
    """xs - An array of frequencies, one dimension"""
    """n1s - An array of population size among the first group, one dimension"""
    """n2s - An array of population size among the second group, one dimension"""

    ns = n1s + n2s
    num = np.multiply(xs, n1s)
    num = np.divide(num, ns)
    denum = n1s * n2s * xs * (1 - xs / ns) / (ns * (ns - 1))
    denum = denum.sum()
    
    if denum == 0:
        return 0
    
    aMin = np.maximum(0, xs - n2s).sum()
    numMin = aMin - num.sum()
    numMin = numMin**2
    aMax = np.minimum(xs, n1s).sum()
    numMax = aMax - num.sum()
    numMax = numMax**2
    
    num = max(numMin, numMax)
    return num / denum

def envelope(xs, n1s, n2s):
    """Given the frequencies in the various populations, compute the envelope"""
    """xs - An array of frequencies, one dimension"""
    """n1s - An array of population size among the first group, one dimension"""
    """n2s - An array of population size among the second group, one dimension"""

    # Otherwise we can compute the enveloppe
    ns = n1s + n2s
    betaLs = np.vstack([xs, n2s]).max(axis = 0)
    betaRs = np.vstack([xs, n1s]).max(axis = 0)
    betaLs = n1s * betaLs /(ns ** 2)
    betaRs = n2s * betaRs /(ns ** 2)
    piL = np.argsort(betaLs)
    piR = np.argsort(betaRs)

    xStars = ns * 1
    pMinLs = np.ones(len(piL,))
    for i in range(len(piL)):
        xStars[piL[i]] = max(xs[piL[i]], n2s[piL[i]])
        pMinLs[i] = minimal_p_value(xStars, n1s, n2s)

    xStars = ns * 1
    pMinRs = np.ones(len(piR,))
    for i in range(len(piR)):
        xStars[piR[i]] = max(xs[piL[i]], n1s[piL[i]])
        pMinRs[i] = minimal_p_value(xStars, n1s, n2s)

    env = np.append(pMinLs, pMinRs).max()
    return env

def old_envelope(xs, n1s, n2s):
    """Given the frequencies in the various populations, compute the envelope"""
    """xs - An array of frequencies, one dimension"""
    """n1s - An array of population size among the first group, one dimension"""
    """n2s - An array of population size among the second group, one dimension"""

    thresh = xs < np.maximum(n1s, n2s)
    # If were are below the threshold, we do not compute the enveloppe
    # since there will be no prunning anyway
    if any(thresh):
        return minimal_p_value(xs, n1s, n2s)
    
    # Otherwise we can compute the enveloppe
    ns = n1s + n2s
    betaLs = n1s * xs /(ns ** 2)
    betaRs = n2s * xs /(ns ** 2)
    piL = np.argsort(betaLs)
    piR = np.argsort(betaRs)

    xStars = ns * 1
    pMinLs = np.ones(len(piL,))
    for i in range(len(piL)):
        xStars[piL[i]] = xs[piL[i]]
        pMinLs[i] = minimal_p_value(xStars, n1s, n2s)

    xStars = ns * 1
    pMinRs = np.ones(len(piR,))
    for i in range(len(piR)):
        xStars[piR[i]] = xs[piR[i]]
        pMinRs[i] = minimal_p_value(xStars, n1s, n2s)

    env = np.append(pMinLs, pMinRs).max()
    return env
