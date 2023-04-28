import numpy as np
from qiskit import QuantumRegister as QReg
from qiskit import QuantumCircuit as QCirc
from qiskit.circuit.library import AND, OR

class QCircH(QCirc):
    '''Helper class for defining quantum circuits'''
    def __init__(self, qc):
        super().__init__(*qc.qregs)
        self.compose(qc, qubits=qc.qubits, inplace=True)

    def find_reg(self, name):
        '''Find the named quantum register'''
        return [reg for reg in self.qregs if reg.name == name][0]

class Paths(QCircH):
    '''When viewed as a Boolean function returns paths given edges
        where paths has True at ix.path(d, j, k) iff 
        there is a path of length <= d+2 between nodes j and k'''
    def __init__(self, ix, num_d_steps):        
        edges = QReg(ix.num_edges, 'edges')
        paths = QReg(num_d_steps*ix.num_paths_per_d, 'paths')
        ancs  = QReg(ix.num_nodes-2, 'ancs')
        qc = QCirc(edges, paths, ancs)

        for d in range(num_d_steps):
            for i in range(ix.num_paths_per_d):
                pairs, single, out = paths_test1_vars(
                    edges, paths, ix, d, i)          
                triples = list(zip(pairs, ancs))
                for (control1, control2), anc in triples: 
                    qc.ccx(control1, control2, anc) # compute and
                qc.append(OR(len(triples)+1).to_gate(), 
                    ancs[:] + [single, out])
                for (control1, control2), anc in triples: 
                    qc.ccx(control1, control2, anc) # uncompute and

        super().__init__(qc)

def paths_test1_vars(edges, paths, ix, d, i):
    '''Setup variables to perform Test 1 using ANDs and ORs'''
    j, k = ix.to_edge(i)
    def prev(p):
        if d == 0:
            return edges[ix.edge(j, p)]
        else:
            return paths[ix.path(d-1, j, p)]
    def next(p):
        return edges[ix.edge(p, k)]
    pairs = [(prev(p), next(p)) for p in ix.nodes_complement([j,k])]
    single = prev(k)
    out = paths[ix.path(d, j, k)]
    return pairs, single, out

class Dist_check(QCircH):
    '''When viewed as a Boolean function returns out given edges
        where out is True iff the graph has the given distances.
        Here ds = [(i0, d0), (i1, d1), ...] and d0 is the distance 
        between nodes j and k where i0 = ix.edge(j, k).'''
    def __init__(self, ix, ds):
        flags = dist_check_test2_flags(ix, ds)
        num_d_steps = flags.shape[0] - 1
        flags = flags.flatten().tolist()

        edges = QReg(ix.num_edges, 'edges')
        paths = QReg(num_d_steps*ix.num_paths_per_d, 'paths')
        ancs  = QReg(ix.num_nodes-2, 'ancs')
        out   = QReg(1, 'out')
        qc = QCirc(edges, paths, ancs, out)

        qc.append(Paths(ix, num_d_steps).to_gate(), 
            edges[:] + paths[:] + ancs[:]) # compute paths
        vals = edges[:ix.num_paths_per_d] + paths[:]
        qc.append(AND(len(vals), flags=flags).to_gate(),
            vals + out[:])
        qc.append(Paths(ix, num_d_steps).to_gate().inverse(), 
            edges[:] + paths[:] + ancs[:]) # uncompute paths

        super().__init__(qc)

def dist_check_test2_flags(ix, ds):
    '''Compute the flags needed to perform Test 2 
        using ANDs and NOTs'''
    d_max = max([d for _, d in ds])
    flags = np.full((d_max, ix.num_paths_per_d), 0)
    for i, d in ds:
        if d == 1:
            flags[0, i] = 1
        else:
            flags[0, i] = -1
            for l in range(d-2):
                flags[1+l, i] = -1
            flags[d-1, i] = 1
    return flags

class Dist_check_groups(QCircH):
    '''When viewed as a Boolean function returns out given edges
        where out is True if the graph has the given distances
        
        Distances are grouped like [(ix0, ds0), (ix1, ds1), ...],
        where ix0 is PathIndexing, ds0 = [(i0, d0), (i1, d1), ...],
        and d0 is the distance between nodes j and k 
        where i0 = ix0.edge(j, k).

        The grouping can be used to control the trade-off between 
        the number of qubits and the size of the circuit.
        '''
    def __init__(self, ds_groups):
        paths_size = max([d*ix.num_paths_per_d 
            for ix, ds in ds_groups for _, d in ds])
        ix, _ = ds_groups[0]

        edges = QReg(ix.num_edges, 'edges')
        paths = QReg(paths_size, 'paths')
        ancs1 = QReg(ix.num_nodes-2, 'ancs1')
        ancs2 = QReg(len(ds_groups), 'ancs2')
        out   = QReg(1, 'out')
        qc = QCirc(edges, paths, ancs1, ancs2, out)

        triples = list(zip(ds_groups, ancs2))
        for (ix, ds), anc2 in triples: 
            dist_check = Dist_check(ix, ds)
            paths_size = dist_check.find_reg('paths').size
            qc.append(dist_check.to_gate(),
                edges[ix.edge_permutation] + paths[:paths_size] 
                + ancs1[:] + [anc2]) # compute dist_check
        qc.append(AND(len(ds_groups)).to_gate(),
            ancs2[:] + out[:])
        for (ix, ds), anc2 in reversed(triples): 
            dist_check = Dist_check(ix, ds)
            paths_size = dist_check.find_reg('paths').size
            qc.append(dist_check.to_gate().inverse(),
                edges[ix.edge_permutation] + paths[:paths_size] 
                + ancs1[:] + [anc2]) # uncompute dist_check
        
        super().__init__(qc)

class Bitflip2Phaseflip(QCircH):
    '''Turn a bit flip oracle into a phase flip oracle'''
    def __init__(self, bitflip):
        qc = QCirc(bitflip.num_qubits)
        out = bitflip.num_qubits - 1 
        qc.x(out)
        qc.h(out)
        qc.append(bitflip.to_gate(), 
            list(range(bitflip.num_qubits)))
        qc.h(out)
        qc.x(out)              
        super().__init__(qc)
