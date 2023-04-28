import numpy as np

class EdgeIndexing():
    def __init__(self, num_nodes, perm=None):
        '''Indexing for edges of a graph with given number of nodes
            Nodes can be rearranged using the given permutation.
            When perm is an int, construct a permutation that 
            shifts indices by the given amount.''' 
        self.num_nodes = num_nodes
        if perm is None:
            perm = 0
        if isinstance(perm, int):
            self.perm = list(range(num_nodes))
            self.perm = self.perm[perm:] + self.perm[:perm]
        else:
            self.perm = perm

    def permute(self, i):
        '''Permuted index of given node'''
        return self.perm.index(i)

    def unpermute(self, i):
        '''Inverse of permute'''
        return self.perm[i]

    @property
    def edge_permutation(self):
        '''Permutation of the edges'''
        ix = EdgeIndexing(self.num_nodes)
        return [ix.edge(*self.to_edge(i)) 
            for i in range(self.num_edges)]

    def num_edges_below_level(self, l):
        '''Number of elements in the strict upper triangle of 
            the incidence matrix above given row'''
        out = 0
        for k in range(l):
            out += self.num_nodes - 1 - k
        return out

    def edge(self, i, j):
        '''Index of the edge between nodes i and j'''
        k, l = self.permute(i), self.permute(j)
        if l < k:
            k, l = l, k
        return self.num_edges_below_level(k) + l - k - 1

    def to_edge(self, ind):
        '''Inverse of edge'''
        for i in range(self.num_nodes-1):
            if ind < self.num_edges_below_level(i+1):
                j = ind - self.num_edges_below_level(i) + i + 1
                return self.unpermute(i), self.unpermute(j)
        raise IndexError('ind out of bounds')

    @property
    def num_edges(self):
        '''Total number of edges'''
        return self.num_edges_below_level(self.num_nodes-1)

    def nodes_complement(self, nodes):
        '''Return the nodes not in the given list of nodes'''
        all_nodes = set(range(self.num_nodes))
        return list(all_nodes.difference(set(nodes)))
    
class PathIndexing(EdgeIndexing):
    def __init__(self, num_nodes, perm=None, num_levels=None): 
        '''Indexing for paths of a graph with given number of nodes
            Here perm is a permutation as in EdgeIndexing, 
            and num_levels is the number of nodes from which 
            the paths originate. This number should be < num_nodes,
            and the default value num_nodes - 1 corresponds to 
            the indexing of all paths. The paths originate from 
            the first nodes listed in perm.'''
        if num_levels is None:
            num_levels = num_nodes - 1
        self.num_levels = num_levels
        super().__init__(num_nodes, perm=perm)

    @property
    def num_paths_per_d(self):
        '''Number of paths per d where d is as in path below'''
        return self.num_edges_below_level(self.num_levels)
    
    def path(self, d, i, j):
        '''Index of the dth path between nodes i and j'''
        return d*self.num_paths_per_d + self.edge(i, j)
