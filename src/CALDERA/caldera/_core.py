################################################################################
# This is the main function where we actually take the pieces from all the
# scripts and put it all together
################################################################################

import pandas as pd
import numpy as np
import caldera.Exploration._ExploreDFS as dfs
import caldera.Exploration._ExploreBFS as bfs
from caldera.Statistics._Pop import *
import pickle
from joblib import Parallel, delayed, parallel_backend
import caldera.Statistics._Helper_Tarone as Helper_Tarone
from caldera.Exploration._Helper_Explo import structure
from caldera.Post._dbgwas_format import create_files

def create_G(loc, C, output, comFile, verbose, restart):
    """
    Create the G structure based on the files from step 1 that were generated
    while building the compacted De Bruijn Graph
    
    Args:
        loc: location of the step 1 folder
        C: number of communities to use for k-means if no comFile
        comFile: an optional file for the community file if computed
        output: where to store the population vector
        verbose: whether to print results
        restart: Whether this is a restart.
        
    Return:
        An object of class structure containing the nodes, the nodes's
        neighbors, the nodes' patterns and the nodes' lengths, as well as the
        samples' community assignements and phenotypes.
    
    Example: 
    """
    if restart:
        G = pickle.load(open(loc + 'save_int_G.obj', 'rb'))
        return G
    
    U = np.loadtxt(loc + "/step1/major_pattern.all_rows.binary", dtype = np.int)
    lengths = np.empty(U.shape[0])
    with open(loc + "/step1/graph.nodes") as fp:
        for cnt, line in enumerate(fp):
            lengths[cnt] = len(str(line)) - len(str(cnt)) - 2
    lengths = list(map(int, lengths))
    lengths = np.array(lengths)
    
    seed = 7361
    cutoff = .9
    if verbose:
        print("Computing the population structure")
    
    if comFile == "None":
        pop = assignPops(U.transpose(), C, cutoff)
        np.savetxt(output + "/pop.txt", pop, fmt = '%d')
    else:
        pop = np.loadtxt(comFile, dtype = np.int)
    # Save pop for next steps.
    pattern = (U == True)
    Pheno = np.loadtxt(loc + "/step1/bugwas_input.id_phenotype",
                       skiprows = 1, usecols = 1, dtype = np.int)
    
    if verbose:
        print("Computing the neighbours")
    edges = pd.read_table(loc + "/step1/graph.edges.dbg", header = None,
                          names = ['From', 'To'], usecols = [0,1])
    edges['min'] = edges.min(axis = 1)
    edges['max'] = edges.max(axis = 1)
    edges.drop(["From", "To"], axis = 1, inplace = True)
    edges.drop_duplicates(inplace = True)
    edges = edges.to_numpy()
    
    neighbours = {k: [] for k in range(np.shape(pattern)[0])}
    for i in range(np.shape(edges)[0]):
        edge = edges[i, :]
        neighbours[int(edge[0])].append(int(edge[1]))
        neighbours[int(edge[1])].append(int(edge[0]))
    
    G = bfs.structure(pop, neighbours, pattern, Pheno, lengths)

    if verbose:
        print("Finished computing the structure", flush = True)
        print("", flush = True)
    
    return G

def CALDERA(loc, output, threads, communities, comFile, Lmax, kmax, sMax = 5, batch_size = 10 ** 10,
            alpha = 10 ** (-8), verbose = True, BFS = True, save_int = False, restart = False):
    """

    Run the Tarone exploration procedure on the graph and returns a set of testable subgraphs.
    
    Args:
        loc: location of the step 1 folder. See details below
        output: Folder where to store the output.
        threads: How many threads to use. Passed to joblib.Parallel
        communities: Number of communities to find. This uses kMeans on the top PCs.
        comFile: an optional file where community assignment is specified. 
                 Either communities or comFile must be specified.
        Lmax: Maximum size of each subgraph in bp.
        sMax: maximum number of stages. Default to 5.
        batch_size: Number of batches to use, to limit the breath of the exploration.
                    Results in hybrid exploration. Not used if BFS = False. 
                    Default to 10 ** 10.
        alpha: value used to control the FWER. Default to 10 ** (-8)
        verbose: whether to be verbose or not. Default to True.
        ## For multiple runs
        save_int: Whether to save the results for a future run
        restart: Whether this is a restart with a higher value of sMax. 
        ## Internal parameters
        BFS: Whether to explore the graph in BFS or DFS. Default to True.
        kmax: maximum value for k to avoid exploring all the graph in case
              everything is significant. In that case increase alpha. 
       
    Returns:
        Nothing. Save the vector of testable subgraphs and their envelopes to
        file. 
    """
    G = create_G(loc, communities, output, comFile, verbose, restart)
    if save_int:
        pickle.dump(G, open(loc + 'save_int_G.obj', 'wb'))
        pickle.dump([sMax, Lmax, alpha], open(loc + 'save_int_param.obj', 'wb'))
                              
    # We compute the set of testable subgraphs based on the envelope
    if BFS:
        R, k0, TH = bfs.explore(G = G, alpha = alpha, kmax = kmax,
                                Lmax = Lmax, verbose = verbose,
                                threads = threads, sMax = sMax,
                                batch_size = batch_size, 
                                restart = restart,
                                save_int = save_int,
                                loc = loc)
    else:
        R, k0, TH = dfs.explore(G = G, alpha = alpha, kmax = kmax,
                                Lmax = Lmax, verbose = verbose,
                                threads = threads, sMax = sMax)

    # We further prune this group based on the minimal p-value
    if verbose:
        print("Actually testing", flush = True)
    if verbose:
        print("Computing the p-values", flush = True)
    pvals = np.array([S.TH(G.Pop, G.Pheno) for S in R])
    # pvals = pvals[:,0]
    if verbose:
        print("Filtering based on k0", flush = True)
    sig_graphs = pvals >= TH
    
    if np.sum(sig_graphs) == 0:
        print("We found no significant subgraphs", flush = True)
        return
    else:
        sigs = R[sig_graphs]
        sigs_pvals = pvals[sig_graphs]

        # Save all the output
        # subgraph_file = open(output + '/Sig_subgraphs.obj', 'wb')
        # pickle.dump(sigs, subgraph_file)
        # subgraph_file.close()
        create_files(loc, output, sigs, sigs_pvals, G, verbose)
    return
