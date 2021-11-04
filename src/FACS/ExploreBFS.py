################################################################################
# Implementation of our algorithm to explore all connected closed
# subgraphs in breadth
################################################################################

import sys, os, tempfile, pickle, shutil
import numpy as np
from scipy.stats import chi2
from joblib import Parallel, delayed, parallel_backend, dump, load
from Envelope import *
from Helper_Explo import *
from Helper_Tarone import *
from Subgraphs import *
    
################################################################################
# This is the part where we explore the graph per say
################################################################################

def explore(G, alpha, kmax, Lmax = 10 ** 7, verbose = False, threads = 1,
            backend = "multiprocessing", old_Env = False):
    """Compute the set R by incremental values of k:
    
    Args:
        G: the graph we want to explore.
        alpha: between 0 and 1. FWER we want to control.
        n1s: An vector of group sizes among the first phenotype.
        n2s: An vector of group sizes among the second phenotype.
        kmax: an integer, the maximum value of k we allow.
        Lmax: maximum size of the subgraph
        verbose: a boolean, default to false.
        threads: Number of cores used in the process
        backend: the backend used for parallelization
    
    Return:
        The array R of testable subgraphs.
    """
    # We avoid computing this too often
    n1s, n2s = G.ns()
    n = int(n1s.sum() + n2s.sum())
    nNodes = G.lengths.shape[0]
    # We store the patterns in a memory-array format
    Patterns = G.pattern
    temp_folder = tempfile.mkdtemp()
    if verbose:
        print("Saving intermediary files in " + temp_folder, flush = True)
        
    filename = os.path.join(temp_folder, 'Patterns.mmap')
    if os.path.exists(filename):
        os.unlink(filename)
    
    _ = dump(Patterns, filename)
    Patterns = load(filename, mmap_mode = 'r+')
    
    Neighbours = G.neighbours
    filename = os.path.join(temp_folder, 'Neighbours.mmap')
    if os.path.exists(filename):
        os.unlink(filename)
    
    _ = dump(Neighbours, filename)
    Neighbours = load(filename, mmap_mode = 'r+')
    
    Lengths = G.lengths
    filename = os.path.join(temp_folder, 'Lengths.mmap')
    if os.path.exists(filename):
        os.unlink(filename)
    
    _ = dump(Lengths, filename)
    Lengths = load(filename, mmap_mode = 'r+')
    
    k0, TH = find_ko(np.array([]), alpha , 1)
    
    # Initiate the list of candidates. We create the closure of every subgraph
    # of size 1
    with parallel_backend(backend, n_jobs=threads):
        chunks = np.array(range(nNodes))
        chunks = np.array_split(chunks, threads)
        C = Parallel(verbose = False)(delayed(start)(
                     chunk, G.Pop, Patterns, Neighbours, Lengths, TH, n1s, n2s, old_Env)
                     for chunk in chunks)
        C = np.concatenate(C)
        
        Env = np.array([S.Env for S in C], dtype = np.float64)
        Pvals = np.array([S.pval for S in C], dtype = np.float64)
        Not_Too_Large = np.array([S.length < Lmax for S in C], dtype = np.bool_)
        
        # We then find the updated k0 we considering only subgraphs of size 1
        k0, TH = find_ko(Pvals, alpha)
        # We clean the list and we save testable subgraphs
        are_Testable = Pvals >= TH
        R = C[are_Testable]
        R_Pvals = Pvals[are_Testable]
        Keep = (Env >= TH) + (Not_Too_Large)
        print('We continue with ' + str(C.shape), flush = True)
        C = C[Keep]
        
        # We then explore all subgraphs of size 2 and so on and so forth
        while C.shape[0] > 0:
            if k0 > kmax:
                if verbose:
                    print('Reached kmax value')
                break
            # Compute all the childrens

            chunks = np.array_split(C, threads)
            Childrens = Parallel(verbose = verbose)(delayed(expand)(
                               chunk, G.Pop, Patterns, Neighbours, Lengths,
                               TH, n1s, n2s, old_Env) for chunk in chunks)
            C = np.concatenate(Childrens)
            if C.shape[0] == 0:
                break

            Env = np.array([S.Env for S in C], dtype = np.float64)
            Pvals = np.array([S.pval for S in C], dtype = np.float64)
            Not_Too_Large = np.array([S.length < Lmax for S in C],
                                     dtype = np.bool_)
            
            # We then find the updated k0 we considering only subgraphs of size 1
            k0, TH = find_ko(np.concatenate([Pvals, R_Pvals]), alpha, k0)
            # We clean the list and we save testable subgraphs
            are_Testable = Pvals >= TH
            R = np.concatenate([R[R_Pvals >= TH], C[are_Testable]])
            R_Pvals = np.concatenate([R_Pvals[R_Pvals >= TH],
                                      Pvals[are_Testable]])
            Keep = (Env >= TH) + Not_Too_Large
            print('We continue with ' + str(C.shape), flush = True)
            C = C[Keep]
        
    # We have finished and we return the testable graphs
    # We clean the temp folder
    try:
        shutil.rmtree(temp_folder)
    except OSError:
        pass

    return R, R_Pvals, k0, TH
