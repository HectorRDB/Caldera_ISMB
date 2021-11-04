################################################################################
# This is the main function where we actually take the pieces from all the
# scripts and put it all together
################################################################################

import numpy as np
from coin.Exploration._graph_Westbrook import *
from coin.Exploration._graphExplo import *
from coin.Exploration._helper_functions import *
from caldera.Statistics._Pop import *
from caldera.Statistics._Envelope import *
import pickle

def create_G(loc, C, output, comFile, verbose):
    """
    Create the G structure based on the files from step 1 that were generated
    while building the compacted De Bruijn Graph
    
    Args:
        loc: location of the step 1 folder
        C: number of communities to use for k-means if no comFile
        comFile: an optional file for the community file if computed
        independently
        output: where to store the output (for the log and
                the population vector)
        verbose: whether to print results
        
    Return:
        An object of class structure containing the nodes, the nodes's
        neighbors, the nodes' patterns and the nodes' lengths, as well as the
        samples' community assignements and phenotypes.
    """
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
        print("Computing the population structure",
              file = open(output + "log","a"))
    
    if comFile == "None":
        pop = assignPops(U.transpose(), C, cutoff)
        np.savetxt(output + "pop.txt", pop, fmt = '%d')
    else:
        pop = np.loadtxt(comFile, dtype = np.int)
    # Save pop for next steps.
    pattern = U
    Pheno = np.loadtxt(loc + "step1/bugwas_input.id_phenotype",
                       skiprows = 1, usecols = 1)
    
    if verbose:
        print("Computing the neighbours", file = open(output + "log","a"))
    edges = np.loadtxt(loc + "/step1/graph.edges.dbg", usecols = (0, 1))
    
    neighbours = dict()
    for i in range(np.shape(pattern)[0]):
        neighbours[i] = []
    for i in range(np.shape(edges)[0]):
        edge = edges[i, :]
        neighbours[int(edge[0])].append(edge[1])
        neighbours[int(edge[1])].append(edge[0])
                           
    G = structure(pop, neighbours, pattern, Pheno, lengths)
    return G

def COIN(loc, output, communities, Lmax, kmax, comFile, verbose, alpha):
    """
    The core algorithm. Run the Tarone exploration on the graph and returns a
    set of testable subgraphs.
    
    Args:
        location of the step 1 folder
        output: where to store the output
        communities: Number of communities to find
        Lmax: Maximum size of each subgraph in bp
        kmax: maximum value for k (to avoid exploring all the graph)
        comFile: an optional file for the community file if computed
        independently
        verbose: whether to be verbose or not
       
    Returns:
        Nothing. Save the vector of testable subgraphs and their envelopes to
        file. 
    """
    G = create_G(loc, communities, output, comFile, verbose)
    n1s, n2s = G.ns()
                               
    # We compute the set of testable subgraphs based on the envelope
    R, MinP, k0, TH = explore(G = G, alpha = alpha, n1s = n1s,
                                      n2s = n2s, Lmax = Lmax, kmax = kmax,
                                      verbose = verbose)
    if verbose:
        print("Finished exploring the graph")
        print("We found " + str(len(R)) + "testable subgraphs")
        print("The values of k0  = {k0} and TH = {TH}".format(k0 = k0, TH = TH))
    
    if verbose:
        print("Actually testing")
    pvals = np.array([S.TH(G.Pop, G.Pheno) for S in R])
    pvals = pvals[:,0]
    sig_graphs = pvals >= TH
    print("We found {n} significant subgraphs".format(n = sum(sig_graphs)))
    sigs = R[sig_graphs]
    sig_pvals = pvals[sig_graphs]
    
    # Save all the output
    subgraph_file = open(output + '/Sig_subgraphs.obj', 'wb')
    pickle.dump(sigs, subgraph_file)
    subgraph_file.close()
    
    nodes = [S.nodes for S in sigs]
    f = open(output + '/Sig_subgraphs_nodes.txt','w')
    for node_list in nodes:
        string = ','.join(map(str, node_list))
        f.write(string + '\n')

    f.close()

    f = open(output + '/Sig_subgraphs_pvals','w')
    for pval in sig_pvals:
        f.write(str(pval) + '\n')
    f.close()

    return
