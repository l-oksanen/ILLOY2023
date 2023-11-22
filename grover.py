import numpy as np
import qiskit as qk

def grover(oracle, state, n_steps = -1):
    '''Constuct a quantum circuit that runs n_steps of Grover's 
        algorithm on the qubits given by the indices in state 
        and using the (phase flip) oracle'''

    if n_steps < 0:
        n_steps = grover_n_steps(len(state))

    qc = qk.QuantumCircuit(oracle.num_qubits)
    for qubit in state: # Initialize
        qc.h(qubit)
    for _ in range(n_steps):
        qc.append(oracle, range(oracle.num_qubits))
        # Diffuser steps, for an explanation, see
        # https://learn.qiskit.org/course/ch-algorithms/grovers-algorithm
        # 1. Apply transformation |s> -> |00..0> (H-gates)
        #   where |s> is the uniform superposition
        for qubit in state:
            qc.h(qubit)
        # 2. Apply transformation |00..0> -> |11..1> (X-gates)
        for qubit in state:
            qc.x(qubit)
        # 3. Do multi-controlled-Z gate
        qc.h(state[-1])
        qc.mct(list(state[:-1]), state[-1])  # multi-controlled-toffoli
        qc.h(state[-1])
        # 4. Apply transformation |11..1> -> |00..0>
        for qubit in state:
            qc.x(qubit)
        # 5. Apply transformation |00..0> -> |s>
        for qubit in state:
            qc.h(qubit)
    qc.measure_all()
    return qc    

def grover_n_steps(n_state_qubits):
    '''Optimal number of Grover iterations assuming that the 
        number of solutions is one'''
    theta = np.arcsin(1/np.sqrt(2**n_state_qubits))
    if theta < np.pi/8:
        return int(np.ceil(np.pi/(4*theta) - 1/2))
    elif theta < np.pi/4:
        return 1
    else: 
        return 0