import numpy as np
from circuits import *

def c_paths(edges, ix, num_d_steps):
    '''Return a boolean vector with True at ix.path(d, j, k) iff 
        there is a path of length <= d+2 between nodes j and k'''
    paths = np.full(num_d_steps*ix.num_paths_per_d, False) 
    for d in range(num_d_steps):
        for i in range(ix.num_paths_per_d):
            pairs, single, _ = paths_test1_vars(
                edges, paths, ix, d, i)
            ands = [elem1 and elem2 for (elem1, elem2) in pairs]
            out = single
            for AND in ands:
                out = out or AND
            if out:
                ind = ix.path(d, *ix.to_edge(i))
                paths[ind] = not paths[ind]
    return paths

def c_dist_check(edges, ix, ds):
    '''Check if the graph has the given distances
        Here ds = [(i0, d0), (i1, d1), ...] and d0 is the distance 
        between nodes j and k where i0 = ix.edge(j, k).'''
    out = True
    flags = dist_check_test2_flags(ix, ds)
    num_d_steps = flags.shape[0] - 1
    flags = flags.flatten()
    paths = c_paths(edges, ix, num_d_steps)
    vals = list(edges[range(ix.num_paths_per_d)]) + list(paths)
    for i in range(len(vals)):
        if flags[i] == 1:
            out = out and vals[i]
        elif flags[i] == -1:
            out = out and not vals[i]
    return out

def c_dist_check_groups(edges, ds_groups):
    '''Check if the graph has the given distances        
        Distances are grouped like [(ix0, ds0), (ix1, ds1), ...],
        where ix0 is PathIndexing, ds0 = [(i0, d0), (i1, d1), ...],
        and d0 is the distance between nodes j and k 
        where i0 = ix0.edge(j, k).
        '''
    out = True
    for ix, ds in ds_groups:
        out = out and c_dist_check(
            edges[ix.edge_permutation], ix, ds)
    return out
