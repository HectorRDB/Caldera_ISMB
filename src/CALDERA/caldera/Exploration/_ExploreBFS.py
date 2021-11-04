################################################################################
# Implementation of our algorithm to explore all connected closed
# subgraphs in breadth
################################################################################

import shutil, pickle
import numpy as np
from scipy.stats import chi2
from joblib import Parallel, delayed, parallel_backend
from caldera.Statistics._Envelope import *
from caldera.Statistics._Helper_Tarone import *
from caldera.Exploration._Helper_Explo import *
from caldera.Exploration._Subgraphs import *
from caldera.Exploration._Light_Subgraphs import *
    
################################################################################
# This is the part where we explore the graph per say
################################################################################

class solutions:
    """Used to store the solutions as we build them:
        Slots
        alpha: The FWER control level
        k: the current value of k0
        TH: the corresponding chi-square
        R: the array of testable subgraphs
        Minps: their minimal p-values
        Lmax: maximum allowed size of subgaph
        kmax: maximum allowed value of k
        sMax: maximum number of stages to explore
        """
    #    
    def __init__(self, alpha, nNodes, Lmax, kmax, sMax):
        if alpha is None:
            self.alpha = None
            self.TH = 1
        else:
            self.alpha = alpha
            self.TH = chi2.isf(alpha, 1)
        self.k0 = 1
        self.nNodes = nNodes
        self.Lmax = Lmax
        self.kmax = kmax
        self.sMax = sMax
        self.R = np.array([], dtype = graph)
        self.Minps = np.array([], dtype = np.float64)
    #    
    def build_first_stage(self, threads, G, n1s, n2s, verbose):
            if verbose:
                print("Building first stage")
            chunks = np.array(range(self.nNodes))
            chunks = np.array_split(chunks, threads)
            C = Parallel(verbose = False)(delayed(start)(
                    chunk, G.Pop, G.pattern, G.neighbours, G.lengths, self.TH, n1s, n2s)
                    for chunk in chunks)
            C = np.concatenate(C)
            #
            Env = np.array([S.Env for S in C], dtype = np.float64)
            Minps = np.array([S.minp for S in C], dtype = np.float64)
            Not_Too_Large = np.array([S.length < self.Lmax for S in C], dtype = np.bool_)
            # If alpha is pre-specified we update k0 we considering only subgraphs of size 1
            # Otherwise we first find alpha before updating k0.
            if self.alpha is None:
                self.alpha = find_alpha(C, G.Pop, G.Pheno, top = 10)
                if verbose:
                    message = ('The value of alpha is {alpha}. '.format(alpha = self.alpha))
                    print(message, flush = True)
                    print('', flush = True)
            k0, TH = find_ko(Minps, self.alpha)
            self.k0 = k0
            self.TH = TH
            # We clean the list and we save testable subgraphs
            are_Testable = Minps >= self.TH
            R = C[are_Testable]
            R = np.array([light_graph(S) for S in R])
            self.R = R
            self.Minps = np.array([S.minp for S in R])
            Keep = (Env >= self.TH) * (Not_Too_Large)
            C = C[Keep]
            if verbose:
                print('New k0 is ' + str(self.k0), flush = True)
            #
            return C, 1
    #
    def build_stage(self, stage, C, threads, G, n1s, n2s, verbose):
        chunks = np.array_split(C, threads)
        Childrens = Parallel(verbose = False)(delayed(expand)(
                             chunk, G.Pop, G.pattern, G.neighbours, G.lengths,
                             self.TH, n1s, n2s) for chunk in chunks)
        C = np.concatenate(Childrens)
        if C.shape[0] != 0:
            Env = np.array([S.Env for S in C], dtype = np.float64)
            Minps = np.array([S.minp for S in C], dtype = np.float64)
            Not_Too_Large = np.array([S.length < self.Lmax for S in C],
                                      dtype = np.bool_)
            # We then find the updated k0 we considering only subgraphs of size 1
            k0, TH = find_ko(np.concatenate([Minps, self.Minps]), self.alpha, self.k0)
            self.k0 = k0
            self.TH = TH
            # We clean the list and we save testable subgraphs
            are_Testable = Minps >= self.TH
            new_R = C[are_Testable]
            new_R = [light_graph(S) for S in new_R]
            self.R = np.concatenate([self.R[self.Minps >= TH], new_R])
            self.Minps = np.concatenate([self.Minps[self.Minps >= TH],
                                        [S.minp for S in new_R]])
            Keep = (Env >= TH) * Not_Too_Large
            if verbose:
                print('New k0 is ' + str(self.k0), flush = True)
            C = C[Keep]
        #
        return C
    
    def next_stage(self, stage, C, threads, G, n1s, n2s, verbose, batch_size, message):
        if C.shape[0] > 0 and stage < self.sMax:
            if self.k0 > self.kmax:
                if verbose:
                    print('Reached kmax value', flush = True)
                return
            n_groups = -(-C.shape[0] // batch_size)
            C = np.array_split(C, n_groups)
            Cs = [None] * n_groups
            for group in range(n_groups):
                new_message = message + ", " + str(group + 1) + "/" + str(n_groups)
                if verbose:
                    print(new_message, flush = True)
                C2 = self.build_stage(stage, C[group], threads, G, n1s, n2s, verbose)
                Cs[group] = self.next_stage(stage + 1, C2, threads, G, n1s, n2s, verbose, batch_size, new_message)
            return np.array([S for sublist in Cs for S in sublist])
        if stage == self.sMax:
                return C
        if C.shape[0] == 0:
            return []
        

def explore(G, alpha, kmax, Lmax = 10 ** 7, sMax = 5, verbose = False,
            batch_size = 10 ** 10, threads = 1, backend = "multiprocessing",
            save_int = False, restart = False, loc = './'):
    """Compute the set R by incremental values of k:
    
    Args:
        G: the graph we want to explore.
        alpha: between 0 and 1. FWER we want to control.
        kmax: an integer, the maximum value of k we allow.
        Lmax: maximum size of the subgraph
        sMax: maximum number of stages
        verbose: a boolean, default to false.
        batch_size: the breath size limit
        threads: Number of cores used in the process
        backend: the backend used for parallelization
        ## For multiple runs
        save_int: Whether to save the results for a future run
        restart: Whether this is a restart with a higher value of sMax. 
        loc: where to save the intermediary files if relevant
    
    Return:
        The array R of testable subgraphs.
    """
    # We avoid computing this too often
    n1s, n2s = G.ns()
    n = int(n1s.sum() + n2s.sum())
    nNodes = G.lengths.shape[0]
    # We store the patterns in a memory-array format
    temp_folder = G.offload(verbose = verbose)
    
    # Initiate the list of candidates. We create the closure of every subgraph
    # of size 1

    if restart:
        sols = pickle.load(open(loc + 'save_int_sols.obj', 'rb'))
        C = pickle.load(open(loc + 'save_int_C.obj', 'rb'))
        old_sMax = sols.sMax
        sols.sMax = sMax
    else:
        sols = solutions(alpha, nNodes, Lmax, kmax, sMax)
        if verbose:
            message1 = ('Starting to explore a graph with ' +
                        '{N} nodes.'.format(N = nNodes))
            if alpha is None:
                message2 = 'No alpha value is given. It will be chosen automatically.'
            else:
                message2 = ('The value of alpha is {alpha}. '.format(alpha = alpha) +
                            'k0 is initialized at 1.')
            print(message1, flush = True)
            print(message2, flush = True)
            print(flush = True)
    with parallel_backend(backend, n_jobs=threads):
        if restart:
            message = '1/1'  + ', 1/1' * (old_sMax - 1)
            C = sols.next_stage(old_sMax, C, threads, G, n1s, n2s, verbose, batch_size, message)
        else:
            C, stage = sols.build_first_stage(threads, G, n1s, n2s, verbose)
            # We then explore all subgraphs of size 2 and so on and so forth
            C = sols.next_stage(stage, C, threads, G, n1s, n2s, verbose, batch_size, '1/1')
    
    # We have finished and we return the testable graphs
    # We clean the temp folder
    try:
        shutil.rmtree(temp_folder)
    except OSError:
        pass
    
    if verbose:
        print("Finished exploring the graph", flush = True)
        print("We found " + str(len(sols.R)) + " potentially testable subgraphs", 
              flush = True)

    if save_int:
        pickle.dump(sols, open(loc + 'save_int_sols.obj', 'wb'))
        pickle.dump(C, open(loc + 'save_int_C.obj', 'wb'))

    return sols.R, sols.k0, sols.TH
