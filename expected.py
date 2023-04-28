import numpy as np
from scipy import sparse as spr
import indexing as ind

def imat2edges(ix, imat):
    '''Convert incidence matrix to indexed edges'''
    edges = np.full(ix.num_edges, False)
    for i in range(ix.num_edges):
        edges[i] = imat[ix.to_edge(i)]
    return edges

def dmat_expected(imat):
    '''Compute distances as a matrix using scipy'''
    return spr.csgraph.dijkstra(spr.csr_matrix(imat), 
        directed=False, unweighted=True)

def dmat2ds(dmat):
    '''Convert distance matrix to the form consumed by 
        Dist_check'''
    n, _ = dmat.shape
    ix = ind.PathIndexing(n)
    return ix, [(i, int(dmat[ix.to_edge(i)])) 
        for i in range(ix.num_edges) 
        if dmat[ix.to_edge(i)] != np.inf]

def dmat2ds_groups(dmat):
    '''Convert distance matrix to the form consumed by 
        Dist_check_groups with a separate group for each node'''
    n, _ = dmat.shape
    ds_groups = []
    for o in range(n):
        ix = ind.PathIndexing(n, perm=o, num_levels=1)
        ds = []
        for k in range(o+1, n):
            d = dmat[o, k]
            if d != np.inf:
                ds.append((ix.edge(o, k), int(d)))
        if ds != []:
            ds_groups.append((ix, ds))
    return ds_groups

def bools2str(bs):
    '''Convert a list of booleans to a bit string'''
    return ''.join([str(int(b)) for b in bs])
