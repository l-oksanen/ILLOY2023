from indexing import *
from circuits import *
from expected import *
from classical import *
import qiskit_aer
import qiskit as qk

def find_bits_in_reg(qc, name):
    '''Find the indices of qubits in named quantum register 
        in given quantum circuit'''
    reg = qc.find_reg(name)
    return [qc.find_bit(reg[j]).index for j in range(reg.size)]

def substr(s, ixs):
    '''Pick substring corresponding to given character indexes'''
    return ''.join([s[i] for i in ixs])

def str2bools(s):
    '''Convert a bit string to a list of booleans'''
    return [bool(int(s[i])) for i in range(len(s))]

def run_circuit(circuit, init_regs=[]):
    '''Run quantum circuit with given named registers 
        initialized with given Boolean vectors
        Here init_regs = [(name0, vec0), (name1, vec1), ...]'''
    qc = qk.QuantumCircuit(circuit.num_qubits)
    for name, vec in init_regs:
        bits = find_bits_in_reg(circuit, name)
        pairs = zip(bits, vec)
        for bit, value in pairs:
            if value:
                qc.x(bit)
    qc.append(circuit, range(circuit.num_qubits))
    qc.measure_all()    
    backend = qiskit_aer.AerSimulator(method='statevector')
    compiled_circuit = qk.transpile(qc, backend)
    result = backend.run(compiled_circuit, shots=1).result()
    counts = result.get_counts(compiled_circuit)
    return list(counts.keys())[0][::-1]

def assert_reg_is_zero(outb, qc, name):
    '''Assert that the named register is full of zeros'''
    bits = find_bits_in_reg(qc, name)
    values = substr(outb, bits)
    print(f'Register {name} has {values = } ({bits = })')
    assert values == len(values) * '0'

def q_paths(edges, ix, num_d_steps):
    '''Compute paths using the quantum algorithm'''
    qc = Paths(ix, num_d_steps)
    outb = run_circuit(qc, init_regs=[('edges', edges)])
    assert_reg_is_zero(outb, qc, 'ancs')
    bits = find_bits_in_reg(qc, 'paths')
    return str2bools(substr(outb, bits))

def q_dist_check(edges, ix, ds):
    '''Check distances using the quantum algorithm'''
    qc = Dist_check(ix, ds)
    outb = run_circuit(qc, init_regs=[('edges', edges)])
    print(f'{outb = }')
    assert_reg_is_zero(outb, qc, 'paths')
    assert_reg_is_zero(outb, qc, 'ancs')
    return str2bools(outb)[-1]

def q_dist_check_groups(edges, ds_groups):
    '''Check grouped distances using the quantum algorithm'''
    qc = Dist_check_groups(ds_groups)
    outb = run_circuit(qc, init_regs=[('edges', edges)])
    print(f'{outb = }')
    assert_reg_is_zero(outb, qc, 'paths')
    assert_reg_is_zero(outb, qc, 'ancs1')
    assert_reg_is_zero(outb, qc, 'ancs2')
    return str2bools(outb)[-1]

def bits2imat(ix, b):
    '''Convert bit string to incidence matrix'''
    imat = np.full((ix.num_nodes,ix.num_nodes), False)
    for k in range(ix.num_edges):
        i, j = ix.to_edge(k)
        val = b[k] == '1'
        imat[i, j] = val
        imat[j, i] = val
    return imat    

def paths2dist(ix, edges, paths):
    '''Compute distances to the last node'''
    ds = np.full(ix.num_nodes-1, np.inf)
    for j in range(len(ds)):
        if edges[ix.edge(j,ix.num_nodes-1)]:
            ds[j] = 1
        else:
            for d in range(ix.num_nodes-2):
                if paths[ix.path(d,j,ix.num_nodes-1)]:
                    ds[j] = d+2
                    break
    return ds

def run_test_dists(fun):
    '''Assert that fun computes the correct distances'''
    ns = [3, 4]
    for n in ns:
        ix = PathIndexing(n, num_levels=1, perm=n-1)
        print(f'Distances for {ix.num_edges = }:')
        for i in range(2**ix.num_edges):
            b = format(i, f'0{ix.num_edges}b')
            imat = bits2imat(ix, b)
            edges = imat2edges(ix, imat)
            paths = fun(edges, ix, ix.num_nodes - 2)
            dist = paths2dist(ix, edges, paths)
            dist_exp = dmat_expected(imat)[:-1,-1]
            if np.array_equal(dist, dist_exp):
                print(f'{b} -> {dist} == {dist_exp}')
            else:
                print(f'{b} -> {dist} != {dist_exp} <<<< INCORRECT')
                assert False

def run_test_dist_check_self(fun, use_groups=False):
    '''Assert that fun, performing the distance check, gives 
        True when given distances corresponding to given edges and 
        False when the distances correspond to another graph.'''
    ns = [3, 4]
    for n in ns:
        ix = EdgeIndexing(n)
        print(f'Distance self check with {ix.num_edges = }:')
        for i in range(2**ix.num_edges):
            b = format(i, f'0{ix.num_edges}b')
            imat = bits2imat(ix, b)
            dmat = dmat_expected(imat)
            edges = imat2edges(ix, imat)
            if use_groups:
                ds_groups = dmat2ds_groups(dmat)
                if len(ds_groups) > 0:
                    retval = fun(edges, ds_groups)
                    assert retval            
            else:
                ixp, ds = dmat2ds(dmat)
                if len(ds) > 0:
                    assert fun(edges, ixp, ds)

def run_test_dist_check_tree(fun, imat_tree, use_groups=False):
    '''Assert that fun, performing the distance check, 
        can be used to find a given tree from its distances'''
    if use_groups:
        dists = [dmat2ds_groups(dmat_expected(imat_tree))]
    else:
        dists = dmat2ds(dmat_expected(imat_tree))
    ix, ds = dmat2ds(dmat_expected(imat_tree))
    edges_tree = imat2edges(ix, imat_tree)
    for i in range(2**ix.num_edges):
        b = format(i, f'0{ix.num_edges}b')
        imat = bits2imat(ix, b)
        edges = imat2edges(ix, imat)
        if fun(edges, *dists) == np.array_equal(edges, edges_tree):
            print(f'{b} passes')
        else:
            print(f'{b} FAILS')
            assert False

def run_test_dist_check_trees(fun_dist_check, use_groups=False):
    '''Assert that fun, performing the distance check, 
        correctly finds a couple of trees from their distances'''
    # Test with
    # 0 - 1
    #   / |
    # 3   2
    imat = np.full((4,4), False)
    imat[0, 1] = True
    imat[1, 2] = True
    imat[1, 3] = True
    run_test_dist_check_tree(fun_dist_check, imat, use_groups)

    # Test with
    # 0 - 1 
    # |   |
    # 3   2
    imat = np.full((4,4), False)
    imat[0, 1] = True
    imat[1, 2] = True
    imat[0, 3] = True
    run_test_dist_check_tree(fun_dist_check, imat, use_groups)

    # Test with
    # 0 - 1 
    #     |
    #     2
    imat = np.full((3,3), False)
    imat[0, 1] = True
    imat[1, 2] = True
    run_test_dist_check_tree(fun_dist_check, imat, use_groups)


# Pytest will collect and run the following functions:

def test_c_dists():
    run_test_dists(c_paths)
def test_c_dist_check_self():
    run_test_dist_check_self(c_dist_check)
def test_c_dist_check_groups_self():
    run_test_dist_check_self(c_dist_check_groups, use_groups=True)
def test_c_dist_check_trees():
    run_test_dist_check_trees(c_dist_check)
def test_c_dist_check_groups_trees():
    run_test_dist_check_trees(c_dist_check_groups, use_groups=True)

def test_q_dists():
    run_test_dists(q_paths)
def test_q_dist_check_trees():
    run_test_dist_check_trees(q_dist_check)
def test_q_dist_check_groups_trees():
    run_test_dist_check_trees(q_dist_check_groups, use_groups=True)
