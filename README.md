Implementation of the algorithm introduced in 

> Joonas Ilmavirta, Matti Lassas, Jinpeng Lu, Lauri Oksanen, Lauri Ylinen.
> _Quantum computing algorithms for inverse problems on graphs and an NP-complete inverse problem_.
> <https://arxiv.org/abs/2306.05253>

## Overview of the algorithm

Given an incidence matrix $E(j, k)$, $j,k=1,\dots,n$ for a graph with $n$ nodes and numbers $d(j, k)$ for $j, k \in B \times B$ for some $B \subset \{1,\dots, n\}$, we want to check if the numbers are distances in the graph. 

We begin by computing the shortest paths in the graph. These are encoded by

$$
P(d, j, k) = \begin{cases}
1 & \text{there is a path of length $\le d$ between nodes $j$ and $k$}
\\
0 & \text{otherwise}
\end{cases} 
$$

We set $P(1, j, k) = E(j, k)$, and construct $P(d, j, k)$ recursively by setting $P(d, j, k) = 1$ if the following test passes:

__Test 1.__ $P(d-1, j, k) = 1$ or there is $p$ distinct from $j$ and $k$ such that $P(d-1, j, p) = 1$ and $E(p, k) = 1$.

Method `c_paths` in [classical.py](classical.py) computes $P$ for testing purposes, and the corresponding quantum circuit is `Paths` in [circuits.py](circuits.py). The classical and quantum implementations both use the method `paths_test1_vars` that sets up variables so that the above Test 1 can be performed using ANDs and ORs. 
 
In practice, we will use a minimal number of indices to encode the edges, rather than the incidence matrix.
In `paths_test1_vars`, the variable `paths[i]`  with index `i = ix.path(d, j, k)` corresponds to $P(d + 2, j, k)$. The translation by two comes from base 0 indexing of Python and from the fact that we don't duplicate $E(j,k)$ in `paths`. The duplication is avoided by using the function `prev(p)` to choose a value either from `paths` or from `edges` where the latter corresponds to $E$.

The distances are checked by 

__Test 2.__ 
$P(d, j, k) = 0$ for $d = 1, \dots, d(j, k) - 1$ and $P(d, j, k) = 1$ for $d = d(j, k)$.

This is implemented by `c_dis_check` in [classical.py](classical.py) for testing purposes, and the corresponding quantum circuit is `Dist_check` in [circuits.py](circuits.py). The classical and quantum implementations both use the method `dist_check_test2_flags` that sets up flag variables so that the above Test 2 can be performed using ANDs and NOTs. 

The method `c_dist_check_groups` calls `c_dist_check` several times. With the grouping given by `dmat2ds_groups` in [expected.py](expected.py), its quantum analogue `Dist_check_groups` can be used to implement the algorithm for which the number of qubits scales as $O(n^2)$. Using `Dist_check` directly without grouping yields the complexity $O(n^3)$, but the resulting simplified algorithm can be more efficient for small networks. The simplified algorithm does not use permuted indices, and a reader interested only in the simplified algorithm can ignore the permutations in [indexing.py](indexing.py).

`Bitflip2Phaseflip` is a generic circuit that converts a bit flip oracle to a phase flip oracle. After this conversion, `Dist_check` and `Dist_check_groups` can be used as oracles in Grover's algorithm, implemented in [grover.py](grover.py). Example driver routines running Grover's algorithm are given in [driver.ipynb](driver.ipynb). 

## Unit tests

Tests that attempt to verify correctness of the oracle are implemented in [test_circuits.py](test_circuits.py). The tests are designed to be run using [pytest](https://pytest.org/) framework.

