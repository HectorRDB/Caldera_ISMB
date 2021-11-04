################################################################################
# The functions below are all used to efficiently explore the graph while 
# computing the fewest possible number of envelopes.
################################################################################

import numpy as np
from scipy.stats import chi2
from Subgraphs import *

################################################################################
# Group of functions linked to closure enumeration
################################################################################

def messages(are_Testable, N, Keep, Not_Too_Large, k):
    """Print messages to stdout
    
    """
    # Get the metrics
    not_keep = sum(Keep == False)
    if not_keep == 0:
        fraction_size = 100
        fraction_prune = 100
        fraction_com = 100
    else:
        fraction_size = 100 * round(
                              (not_keep - sum(Not_Too_Large)) / not_keep, 2)
        fraction_prune = 100 * round(
                    sum((Keep == False) * (1 - are_Testable)) / not_keep, 2)
        
    message1 = ('{n1} out of '.format(n1 = len(Keep)) +
                '{n2} subgraphs explored at that '.format(n2 = N) +
                'step are currently potentially testable using the old k0')
    
    message2 = ('{n} subgraphs were pruned, '.format(n = not_keep) +
        'including {size}% because of size, '.format(size = fraction_size) +
        '{prune}% because of the pruning.'.format(prune = fraction_prune))
    
    print(message1, flush = True)
    print('New k0 is ' + str(k), flush = True)
    print(message2, flush = True)
    print(flush = True)
    
    return

def start(nodes, Pop, Patterns, Neighbours, Lengths, TH, n1s, n2s, old_Env):
    """Create the graph from the nodes and compute their envelope in chunks
    
    Args:
        nodes: a list of nodes numbers
        Pop: the population assignement of the samples
        Patterns: The nodes' patterns
        Neighbours: The nodes' neighbours
        Lengths: The nodes' lengths
        TH: the initial value of the threshold
        n1s: An vector of group sizes among the first phenotype.
        n2s: An vector of group sizes among the second phenotype.
    Return:
        The graphs, their envelopes and some info for pruning
    """
    R = [create_Graph(s, Pop, Patterns, Neighbours, Lengths, TH, n1s, n2s, old_Env)
         for s in nodes]
    R = np.array([S for S in R if S.closed and not S.to_prune], dtype = graph)
    return R

def create_Graph(s, Pop, Patterns, Neighbours, Lengths, TH, n1s, n2s, old_Env):
    """Create and return a (closed) subgraph starting from node s of graph G

    Args:
        s: The node
        Pop: the population assignement of the samples
        Patterns: The nodes' patterns
        Neighbours: The nodes' neighbours
        Lengths: The nodes' lengths
        TH: the initial value of the threshold
        n1s: An vector of group sizes among the first phenotype.
        n2s: An vector of group sizes among the second phenotype.
    Return:
        The closure
    """
    S = graph()
    S.new_graph(s, Pop, Patterns, Neighbours, Lengths, TH, n1s, n2s, old_Env)
    return S

def expand(Ss, Pop, Patterns, Neighbours, Lengths, TH, n1s, n2s, old_Env):
    """Find all parents

    Args:
        Ss: an array of objects of class graph for which we want to compute the
        envelope
        Pop: the population assignement of the samples
        Patterns: The nodes' patterns
        Neighbours: The nodes' neighbours
        Lengths: The nodes' lengths
        TH: the initial value of the threshold
        n1s: An vector of group sizes among the first phenotype.
        n2s: An vector of group sizes among the second phenotype.
        
    Return:
        An array of childrens
    """
    if Ss.shape[0] == 0:
        return np.array([], dtype = graph)
    else:
        Childrens =  np.concatenate([S.Childrens(Pop, Patterns, Neighbours,
                                               Lengths, TH, n1s, n2s, old_Env) for
                                   S in Ss])
        return Childrens

################################################################################
# The structure stores all the info we need on the graph: the list of nodes, the
# list of neighbours (computed from the list of edges), the vector of
# populations and of phenotype, and the lengths of each node.
################################################################################

class structure:
    """Contains the main initial files for the analysis:
        Slots
        Pop: The vector of population
        neighbors: the list of edges
        pattern: the matrix of presence / absence of nodes by samples
        Pheno: the binary vector of categories
        lengths: the size in bp of the nodes
        """
    def __init__(self, Pop, neighbours, nodes, pheno, lengths):
        self.Pop = Pop
        self.neighbours = neighbours
        self.pattern = nodes
        self.Pheno = pheno
        self.lengths = lengths
    
    def ns(self):
        # We have k groups
        # We returns two vectors of length k with the number of sample in each
        # group, split by phenotype.
        clusters = np.unique(self.Pop)
        clusters.sort()
        n1s = np.zeros(len(clusters))
        n2s = np.zeros(len(clusters))
        
        for clus in clusters:
            n1s[clus] = sum(self.Pheno[self.Pop == clus])
            n2s[clus] = sum(self.Pop == clus) - n1s[clus]
        
        return n1s, n2s
