import numpy as np
import pandas as pd
import copy
from caldera.Statistics._Envelope import *


class graph:
    """A graph is a subgraph of the De Bruijn graph, represented in a block-forest structure.
    
    Slots:
        ys: the presence/absence pattern of the subgraph among the samples.
        neighbours: the set of nodes that are the neighbours of the subgraph in the DBG.
        length: the number of bp in the subgraph
        nodes: the list of current nodes.
        Env: the subgraph's envelope
        minp: the subgraph's minimal p-value
        max: the maximum node value
        closed: whether the subgraph is closed or not
    """
    
    def __init__(self):
        """Creating the closure of the subgraph {node}
        
        Return:
            An empty subgraph.
        """
        # Create the subgraph with one node
        ## List of presence / absence in the subgraph for each sample
        self.ys = np.empty((1,))
        self.common = np.empty((1,))
        self.diff_items = []
        self.diff_patterns = np.empty((1,))
        self.iS = -1
        self.neighbours = []
        ## List of neighoring nodes that can be added
        ## Length in bp of all the nodes
        self.length = 0
        self.nodes = np.array([])
        self.Env = 0
        self.minp = 0
        self.prunable = False
        self.to_prune = False
        self.max = 0
        self.closed = True
        return
    
    def new_graph(self, node, Pop, Patterns, Neighbours, Lengths, TH, n1s, n2s):
        """Creating the closure of the subgraph {node}
        
        Args:
            node: The node from which to start
            Pop: the population assignement of the samples
            Patterns: The nodes' patterns
            Neighbours: The nodes' neighbours
            Lengths: The nodes' lengths
            TH: the initial value of the threshold
            n1s: An vector of group sizes among the first phenotype.
            n2s: An vector of group sizes among the second phenotype.
        Return:
            the subgraph
        """
        # Create the subgraph with one node
        ## List of presence / absence in the subgraph for each sample
        self.ys = np.asarray(Patterns[node])
        self.compute_envelope(Pop, n1s, n2s)
        # If we can prune directly, we prune
        if self.Env < TH:
            self.to_prune = True
            return
        ## List of neighoring nodes that can be added
        self.neighbours = copy.copy(Neighbours[node])
        ## Length in bp of all the nodes
        self.length = np.asarray(Lengths[node])
        self.nodes = {node}
        self.max = node
        # Close and get common items
        self.close_init(Patterns, Neighbours, Lengths)
        self.common = np.all(Patterns[np.array(list(self.nodes)),], axis = 0)
        self.parent_neighbours = self.neighbours
        if not all(self.ys == self.common):
            self.closed = False
        return
            
    def close_init(self, Patterns, Neighbours, Lengths, L_Stop = 10 ** 5):
        """Creating the closure of a subgraph
           
        Args:
           self: the current subgraph
           Patterns: The nodes' patterns
           Neighbours: The nodes' neighbours
           Lengths: The nodes' lengths

        Return:
           the closure of the subgraph
        """
        closure = (np.asarray(Patterns[self.neighbours,]) | self.ys) == self.ys
        closure = closure.all(axis = -1 )
        closure = np.where(closure)[0]
        if closure.shape[0] == 0:
            return
        closed_neighbours = [ne for i, ne in enumerate(self.neighbours) if i in closure]
        while len(closed_neighbours) > 0:
            ne = closed_neighbours.pop()
            if ne > self.max:
                self.closed = False
                break
            new_closed_neighbours = self.add_closure_node(
                                            ne, Patterns, Neighbours, Lengths)
            closed_neighbours.extend(new_closed_neighbours)
            if self.length > L_Stop:
                self.closed = (len(closed_neighbours) == 0)
                break
        return
       
    def add_closure_node(self, node, Patterns, Neighbours, Lengths):
        """ Add a node to the growing subgraph that respect the closure

        Args:
           self: the current subgraph
           node: the node to add
           Patterns: The nodes' patterns
           Neighbours: The nodes' neighbours
           Lengths: The nodes' lengths

        Return:
           S with the new node.
        """
        if node not in self.neighbours:
            raise NameError('This node should not be added since it is not a neighbour')
        self.neighbours.remove(node)
        self.nodes.add(node)
        self.length = self.length + Lengths[node]
        # We add only the new neighbours
        new_neighbours = [ne for ne in Neighbours[node] if ne not in
                         self.nodes and ne not in self.neighbours]
        self.neighbours.extend(new_neighbours)
        new_closure = (np.asarray(Patterns[new_neighbours,]) | self.ys) == self.ys
        new_closure = np.all(new_closure, axis = 1)
        new_closure = np.where(new_closure)[0]
        new_closed_neighbours = [ne for i, ne in enumerate(new_neighbours) if
                                 i in new_closure]
        return new_closed_neighbours
            
    def Childrens(self, Pop, Patterns, Neighbours, Lengths, TH, n1s, n2s,
                  L_Stop = 10 ** 5):
        """ Return the closure of all self + v for all neighbours v of self
        
        Args:
            self: the current subgraph
            Pop: the population assignement of the samples
            Patterns: The nodes' patterns
            Neighbours: The nodes' neighbours
            Lengths: The nodes' lengths
            TH: the initial value of the threshold
            n1s: An vector of group sizes among the first phenotype.
            n2s: An vector of group sizes among the second phenotype.
            L_Stop: maximum size of the closure before we stop
        
        Return:
            an array of graph, all children of the subgraph
        """
        # Get the number of neighbours
        if len(self.neighbours) == 0:
            childrens = np.empty(0, dtype = graph)
            return childrens
        closure = np.asarray(Patterns[self.neighbours,]) | self.ys
        closure, indexes = np.unique(closure, axis = 0, return_index = True)
        iSs = np.empty(len(indexes), dtype = graph)
        for i, index in enumerate(indexes):
            ne = self.neighbours[index]
            iSs[i] = compute_iS(self, ne, Patterns)[0]
        # Compute all direct children and siblings
        children_indexes = [index for i, index in enumerate(indexes) if
                            iSs[i] > self.iS and not self.ys[iSs[i]]]
        children_closures = closure[[i for i in range(len(indexes)) if
                            iSs[i] > self.iS and not self.ys[iSs[i]]], ]
        children_iS = iSs[[i for i in range(len(indexes)) if
                           iSs[i] > self.iS and not self.ys[iSs[i]]]]
        childrens = enumerate_childrens(self, children_indexes, children_closures,
                                        children_iS, Patterns, Neighbours,
                                        Lengths, TH, Pop, n1s, n2s, L_Stop)
        siblings_closure = closure[iSs == self.iS, ]
        siblings_indexes = indexes[iSs == self.iS, ]
        siblings = enumerate_siblings(self, siblings_indexes, siblings_closure,
                                      Patterns, Neighbours, Lengths, TH, Pop,
                                      n1s, n2s, L_Stop)
        return np.concatenate([childrens, siblings])
    
    def add_parent_node(self, node, Patterns, common, Neighbours, Lengths, iS,
                        new_diff_patterns):
        """ Add a node to the growing subgraph
         
        Args:
            self: the current subgraph
            node: the node to add
            Patterns: The nodes' patterns
            common: the subgraph common objects after adding the node
            Neighbours: The nodes' neighbours
            Lengths: The nodes' lengths
            iS: the new discriminant item
            new_diff_patterns: the new patterns that we must not enumerate

        Return:
            The subgraph with the new node.
        """
        if node not in self.neighbours:
            raise NameError('This node should not be added since it is not a neighbour')
        # We add the node
        self.max = max(self.nodes)
        self.iS = iS
        self.ys = self.ys | np.asarray(Patterns[node,])
        self.common = common
        self.neighbours.remove(node)
        self.nodes.add(node)
        self.length = self.length + Lengths[node]
        # We remember the neighbours for condition 3
        neighbours_patterns = np.asarray(Patterns[self.neighbours, iS])
        # We check condition 3
        self.parent_neighbours = [ne for i, ne in enumerate(self.neighbours) if
                                  not neighbours_patterns[i]]
        self.parent_neighbours = set(self.parent_neighbours)
        # We add only the new neighbours
        new_neighbours = [ne for ne in Neighbours[node] if
                          ne not in self.nodes and ne not in self.neighbours]
        self.neighbours.extend(new_neighbours)
        self.diff_items = np.where(self.ys == False)[0]
        self.diff_patterns = new_diff_patterns[:, self.diff_items]
        return
    
    def add_sibling_node(self, node, pattern, common, Neighbours, Lengths,
                         new_diff_patterns, iS_node):
        """ Add a node to the growing subgraph without change iS

        Args:
           self: the current subgraph
           node: the node to add
           pattern: The new pattern of the subgraph
           common: the subgraph common objects after adding the node
           Neighbours: The nodes' neighbours
           Lengths: The nodes' lengths
           new_diff_patterns: the new patterns that we must not enumerate
           iS_node: whether the new node contains iS or not

        Return:
            Whether or not we can add this node.
            The subgraph with the new node.
        """
        if node not in self.neighbours:
            raise NameError('This node should not be added since it is not a neighbour')
        if node in self.parent_neighbours:
            return True
        if not iS_node and node > self.max:
            return True
        self.ys = pattern
        reduced_pattern = pattern[self.diff_items]
        explored = (self.diff_patterns | reduced_pattern) == reduced_pattern
        explored = np.all(explored, axis = -1)
        if np.any(explored):
            return True
        self.common = common
        self.neighbours.remove(node)
        self.nodes.add(node)
        self.length = self.length + Lengths[node]
        # We add only the new neighbours
        new_neighbours = [ne for ne in Neighbours[node] if
                          ne not in self.nodes and ne not in self.neighbours]
        self.neighbours.extend(new_neighbours)
        # We remember the diff_items and diff_patterns to avoid redundancy
        diff_items = np.where(pattern == False)[0]
        diff_positions = np.array([i for i, obj in enumerate(self.diff_items)
                                   if obj in diff_items], dtype = int)
        self.diff_patterns = np.vstack([self.diff_patterns[:,diff_positions],
                                        new_diff_patterns[:,diff_items]])
        self.diff_items = diff_items
        return False

    def close(self, Patterns, Neighbours, Lengths, L_Stop = 10 ** 5):
        """Creating the closure of a subgraph
           
        Args:
           self: the current subgraph
           Patterns: The nodes' patterns
           Neighbours: The nodes' neighbours
           Lengths: The nodes' lengths

        Return:
           the closure of the subgraph
        """
        closure = (np.asarray(Patterns[self.neighbours,]) | self.ys) == self.ys
        closure = np.all(closure, axis = 1)
        closure = np.where(closure)[0]
        if closure.shape[0] == 0:
            return
        closed_neighbours = [ne for i, ne in enumerate(self.neighbours) if i in closure]
        all_closed_neighbours = copy.copy(closed_neighbours)
        ob_items = np.where(self.common)[0]
        ob_items = [item for item in ob_items if item > self.iS]
        while len(closed_neighbours) > 0:
            ne = closed_neighbours.pop()
            # Condition 1
            if len(ob_items) > 0 and not Patterns[ne, ob_items].all():
                self.closed = False
                break
            # Condition 2
            if not Patterns[ne, self.iS] and ne > self.max:
                self.closed = False
                break
            # Condition 3
            if ne in self.parent_neighbours:
                self.closed = False
                break
            new_closed_neighbours = self.add_closure_node(
                                            ne, Patterns, Neighbours, Lengths)
            closed_neighbours.extend(new_closed_neighbours)
            all_closed_neighbours.extend(new_closed_neighbours)
            if self.length > L_Stop:
                self.closed = (len(closed_neighbours) == 0)
                break
        if self.closed and np.any(self.common):
            new_commons = np.all(Patterns[all_closed_neighbours,], axis = 0)
            self.common = self.common & new_commons
        return
   
    
    ############################################################################
    # The functions below (frequencies, TH) are all focused on
    # whether the graph is testable or actually significant. This is the part
    # that will be helpfull for Tarone'strick later on.
    ############################################################################
    
    def frequencies(self, Pop):
        """ Return the margins for the first column of the tables"""
        clusters = np.unique(Pop)
        clusters.sort()
        xs = np.zeros(len(clusters))
        for clus in clusters:
            xs[clus] = self.ys[Pop == clus].sum()
        return xs
    
    def compute_envelope(self, Pop, n1s, n2s):
        """ Compute the envelope of a subgraph"""
        freqS = self.frequencies(Pop)
        thresh = freqS >= np.maximum(n1s, n2s)
        self.prunable = all(thresh)
        self.Env = envelope(freqS, n1s, n2s)
        self.minp = minimal_p_value(freqS, n1s, n2s)
            
    def TH(self, Pop, Pheno):
        """ Compute the p-value fo the subgraph (up to a chi-square transformation)"""
        clusters = np.unique(Pop)
        clusters.sort()
        n1s = np.zeros(len(clusters))
        n2s = np.zeros(len(clusters))
        a_s = np.zeros(len(clusters))
        xs = np.zeros(len(clusters))
        
        for clus in clusters:
            a_s[clus] = self.ys[(Pop==clus) & (Pheno==1)].sum()
            xs[clus] = self.ys[Pop == clus].sum()
            n1s[clus] = sum(Pheno[Pop == clus])
            n2s[clus] = sum(Pop == clus) - n1s[clus]
        
        return Th(a_s, xs, n1s, n2s)

def enumerate_childrens(S, children_indexes, children_closures, children_iS,
                        Patterns, Neighbours, Lengths, TH, Pop, n1s, n2s, L_Stop):
    childrens = np.empty(len(children_indexes), dtype = graph)
    for i, index in enumerate(children_indexes):
        S2 = copy.deepcopy(S)
        ne = S.neighbours[index]
        iS, pattern, common, iS_node = compute_iS(S2, ne, Patterns)
        new_diff_patterns = children_closures[[j for j in
            range(len(children_indexes)) if j > i and children_iS[j] == iS],]
        S2.add_parent_node(ne, Patterns, common, Neighbours, Lengths, iS,
                           new_diff_patterns)
        # If we can prune directly, we prune
        S2.compute_envelope(Pop, n1s, n2s)
        if S2.Env < TH:
            continue
        # We close while checking conditions 2 and 3
        S2.close(Patterns, Neighbours, Lengths, L_Stop)
        childrens[i] = S2
    childrens = childrens[childrens != None]
    childrens = np.array([S2 for S2 in childrens if S2.closed])
    return childrens

def enumerate_siblings(S, siblings_indexes, siblings_closure, Patterns,
                       Neighbours, Lengths, TH, Pop, n1s, n2s, L_Stop):
    siblings = np.empty(len(siblings_indexes), dtype = graph)
    fam_size = siblings.shape[0]
    for i, index in enumerate(siblings_indexes):
        S2 = copy.deepcopy(S)
        ne = S.neighbours[index]
        iS, pattern, common, iS_node = compute_iS(S2, ne, Patterns)
        explored = S2.add_sibling_node(ne, pattern, common, Neighbours, Lengths,
                                       siblings_closure[(i+1):fam_size, ], iS_node)
        if explored:
            continue
        # If we can prune directly, we prune
        S2.compute_envelope(Pop, n1s, n2s)
        if S2.Env < TH:
            continue
        # We close while checking conditions 2 and 3
        S2.close(Patterns, Neighbours, Lengths, L_Stop)
        siblings[i] = S2
    siblings = siblings[siblings != None]
    siblings = np.array([S2 for S2 in siblings if S2.closed])
    return siblings


def compute_iS(S, node, Patterns):
    if node not in S.neighbours:
        raise NameError('This node should not be added since it is not a neighbour')
    pattern = np.asarray(Patterns[node])
    common = S.common & pattern
    pattern = S.ys | pattern
    iS = np.where(pattern & (common == False))[0]
    iS = iS.max()
    return iS, pattern, common, Patterns[node, iS]
