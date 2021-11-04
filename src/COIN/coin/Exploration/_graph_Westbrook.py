################################################################################

# Current implementation of the graph structure described in the paper by
# Jeffery Westbrook and Robert E. Tarjan: Maintaining bridge-connected and
# biconnected components on-line, 1992.
# It may be merged with the graphTools.py file eventually

################################################################################

import numpy as np
import igraph, sys, os
from caldera.Statistics._Envelope import *

################################################################################
# The subgraphs are made of two blocks: circle nodes and square nodes
# The circle nodes are not real nodes, they represent the blocks and are
# represented by negative numbers (and plotted in blue). The square nodes are
# the real nodes from the subgraphs and have postive numbers. They are plotted
# in red.
################################################################################

class circleNode:
    """circleNodes represent blocks.
    
    A label of a circle node is a strictly negative number. We incremement by -1. When we fuse several blocks,
    we take the highest value (i.e the one closer to zero).
    """

    def __init__(self, label, parentLabel = None):
        """Create a circle node:
        
        Args:
            label: a negative integer. The node label
            parentLabel: a positive integer. The parent of a circle node can either be a square node or NULL.
            
        Returns:
            The circleNode
            
        The children of a circle nodes are all square nodes.
        """
        self.label = label
        self.parentLabel= parentLabel
        self.children = []

class squareNode:
    """squareNode represent real vertices. 
    
    They are labelled by positive integers, from the De-Bruijn graph."""
    
    def __init__(self, label, parentLabel):
        """Create a circle node:
        
        Args:
            label: A positive integer. The node label from the DBG.
            parentLabel: A negative interger. A parent of a square node is always a circle node.
        
        Returns:
            The squareNode
            
        The blocks of a square nodes are all circle nodes.
        """
        
        self.label = label
        self.parentLabel= parentLabel
        self.blocks = [parentLabel]

class graph:
    """A graph is a subgraph of the De Bruijn graph, represented in a block-forest structure.
    
    Slots:
        circleNodes: a dictionary of the blocks.
        squareNodes: a dictionary of the real vertices in the subgraph.
        articulations: the subset of squareNodes that are articulation points
        ys: the presence/absence pattern of the subgraph among the samples.
        neighbours: the set of nodes that are the neighbours of the subgraph in the DBG.
        length: the number of bp in the subgraph
        futureLabel: the name of the next block to be created.
    """
    
    def __init__(self, n):
        """Create a subgraph.
        
        Args:
            n: the number of samples
        
        Return:
            An empty subgraph. 
        """
        # Used to define the labels of the circle nodes
        self.futureLabel = -1
        # List of all the "real" nodes
        self.squareNodes = dict()
        # List of all the blocks
        self.circleNodes = dict()
        # Subset of the square nodes
        self.articulations = dict()
        # List of presence / absence in the graph for each sample
        self.ys = np.zeros(n, dtype=bool)
        # List of neighoring nodes that can be added
        self.neighbours = np.array([])
        # Length in bp of all the nodes
        self.length = 0

    ############################################################################
    # The functions below (add_node, condense_path and find_path) are all
    # focused on how to add a node to the subgaph while maintaining block
    # structure. This is where we rely on the work of Westbrook and Trajan .
    ############################################################################
    
    def add_node(self, nodeNumber, pattern, neighbours, length):
        """ Add a node to the growing subgraph
        
        Args:
            nodeNumber: the label of the node to add
            pattern: the presence / absence pattern of this node
            neighbours: the list of neighbours of the node
            length: the length in bp of the node
        
        Return:
            the same subgraph with a new node, andd accordingly updated slots.
        """
        # First we add the node to the vectors of ys
        self.ys = self.ys | (pattern > 0)
        
        # Then we expand the length
        self.length = self.length + length
        
        # Then we create the new node with a circleNode as parent
        blockInd = self.futureLabel
        self.futureLabel = self.futureLabel - 1
        newCircleNode = circleNode(label = blockInd)
        self.circleNodes[blockInd] = newCircleNode
        
        # Add the new node to the graph
        newSquareNode = squareNode(label = nodeNumber,
                                   parentLabel  = blockInd)
        self.squareNodes[nodeNumber] = newSquareNode
        self.circleNodes[blockInd].children = [nodeNumber]
        
        # Add edges between this new node and all its neigbors already in the
        # graph, potentialy merging some blocks along the way. Also add the
        # other neighbours to the list of neighbours.
        for neighbour in neighbours:
            blockInd = self.squareNodes[nodeNumber].parentLabel
            # If the neighbour is in the graph
            if neighbour in self.squareNodes.keys():
                # If the two nodes are already on the same block, nothing to do
                # here. This is case 1 from the figure
                if blockInd in self.squareNodes[neighbour].blocks:
                    continue
                elif (self.circleNodes[blockInd].parentLabel == None and 
                              len(self.circleNodes[blockInd].children) == 1):
                # We are in case 2 from the figure. 
                    if self.squareNodes[neighbour].blocks == [-1]:
                    # If this node we add is only the second node, we reroot the
                    # graph and delete the -1 node
                        self.squareNodes[neighbour].parentLabel = blockInd
                        self.squareNodes[neighbour].blocks = [blockInd]
                        self.circleNodes[blockInd].children.append(neighbour)
                        self.circleNodes.pop(-1, None)
                    else:
                        # Otherwise we name that neighbour the parent of the node 
                        # and add it to the list of articulation points
                        self.circleNodes[blockInd].parentLabel= neighbour
                        self.articulations[neighbour] = self.squareNodes[neighbour]
                        self.squareNodes[neighbour].blocks.append(blockInd)
                else:
                # We are in case 3. We find the path between the two blocks and
                # we merge everything along the way
                    self.condense_path(nodeNumber, neighbour)
            else:
                # We add the node's neighour to the list of neighbour if it is
                # not already there.
                if neighbour not in self.neighbours:
                    self.neighbours = np.append(self.neighbours, neighbour)
        # Remove the node from the list of neighbours
        self.neighbours = self.neighbours[np.where(self.neighbours != nodeNumber)]
    
    def condense_path(self, nodeNumber, neighbour):
        path, root = self.find_path(nodeNumber, neighbour)
        circlePath = [node for node in path if node < 0]
        # If the root is a square node, we merge all circle nodes. All children
        # of those nodes are rerooted to this new node.
        if root >= 0:
            rootCircleNodeLabel = max(circlePath)
            rootCircleNode = circleNode(label = rootCircleNodeLabel,
                                        parentLabel = root)
        else:
            rootCircleNodeLabel = root
            rootCircleNode = self.circleNodes[root]
        
        # Reroot all children to that node
        for node in circlePath:
            # If the node is a square, just reroot it and change its blocks
            for sNode in self.circleNodes[node].children:
                isArticulation = len(self.squareNodes[sNode].blocks) > 1
                self.squareNodes[sNode].parentLabel = rootCircleNodeLabel
                blocks = self.squareNodes[sNode].blocks
                blocks = [block for block in blocks if block not in circlePath]
                self.squareNodes[sNode].blocks = blocks
                self.squareNodes[sNode].blocks.append(rootCircleNodeLabel)
                if not sNode in rootCircleNode.children:
                    rootCircleNode.children.append(sNode)
                # Remove the node from the list of articulations if it is not
                # one anymore
                if isArticulation and len(self.squareNodes[sNode].blocks) == 1:
                    del self.articulations[sNode]
            del self.circleNodes[node]
            # Take care of the root if it is a square node
        if root >= 0:
            blocks = self.squareNodes[root].blocks
            blocks = [block for block in blocks if block not in circlePath or
                                                  block == rootCircleNodeLabel]
            self.squareNodes[root].blocks = blocks
        # Add the new circle node
        self.circleNodes[rootCircleNodeLabel] = rootCircleNode
    
    def find_path(self, leftNode, rightNode):
        # Initiate the search by looking at the parent of each node.
        leftPointer = self.squareNodes[leftNode].parentLabel
        leftPath = [leftNode, leftPointer]
        rightPointer = self.squareNodes[rightNode].parentLabel
        rightPath = [rightNode, rightPointer]
        
        # The path up the tree will alternate between circle and square nodes.
        leftParent = "circle"
        rightParent = "circle"
        
        # We go up the tree until we find a common node.
        while not ((leftPointer in rightPath) | (rightPointer in leftPath)):

            if leftParent == "circle":
                if leftPointer != None:
                    leftPointer = self.circleNodes[leftPointer].parentLabel
                    if leftPointer != None:
                        leftPath.append(leftPointer)
                        leftParent = "square"
            else:
                leftPointer = self.squareNodes[leftPointer].parentLabel
                leftPath.append(leftPointer)
                leftParent = "circle"
            
            if rightParent == "circle":
                if rightPointer != None:
                    rightPointer = self.circleNodes[rightPointer].parentLabel
                    if rightPointer != None:
                        rightPath.append(rightPointer)
                        rightParent = "square"
            else:
                rightPointer = self.squareNodes[rightPointer].parentLabel
                rightPath.append(rightPointer)
                rightParent = "circle"
        
        # We now just need to reconstruct the path
        if leftPointer in rightPath:
            j = rightPath.index(leftPointer)
            path = leftPath + rightPath[0:j]
            newRoot = leftPointer
        else:
            j = leftPath.index(rightPointer)
            path = leftPath[0:j] + rightPath
            newRoot = rightPointer
        return path, newRoot
    
    ############################################################################
    # The functions below (frequencies, TH, candidates) are all focused on
    # whether the graph is testable or actually significant. This is the part
    # that will be helpfull for Tarone'strick later on.
    ############################################################################
    
    def compute_envelope(self, Pop, n1s, n2s):
        """ Compute the envelope of a subgraph"""
        freqS = self.frequencies(Pop)
        thresh = freqS >= np.maximum(n1s, n2s)
        self.prunable = all(thresh)
        self.Env = envelope(freqS, n1s, n2s)
        self.pval = minimal_p_value(freqS, n1s, n2s)
    
    def frequencies(self, Pop):
        """ Return the margins for the first column of the tables"""
        clusters = np.unique(Pop)
        clusters.sort()
        xs = np.zeros(len(clusters))
        for clus in clusters:
            xs[clus] = self.ys[Pop == clus].sum()
        return xs
    
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
    
    def candidates(self, j):
        """ Return whether a node should be added to the graph"""
        # Get a list of all articulation nodes
        articulations = list(self.articulations.keys())
        articulations = list(map(int, articulations))
        
        # Get the list of all nodes in the graph
        nodes = list(self.squareNodes.keys())
        nodes = list(map(int, nodes))
        
        # Get the maximum index of nonVertexSeparator nodes
        nonVertexSeparators = np.array(list(set(nodes).difference(set(articulations))))
        nonVertexSeparators = nonVertexSeparators[nonVertexSeparators != j]
        i = max(nonVertexSeparators, default = -1)
        
        # Get the list of candidates
        return j > i
        
