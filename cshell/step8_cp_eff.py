"""Step 8 - obtain CP-EFF - control-energy reduction.

1:1 port of `[5] Obtain IER_CS.ipynb`.

Random search: starting from driver = CP-IND u CP-EXT u CP-DRV (part1+part2),
sample additional regulator nodes and keep the random subset that yields
the smallest control energy E_u (Yan et al. expression, see paper).
A 100-iteration random search is the default, matching the notebook.

Output file (paper terminology)
-------------------------------
CP-EFF.txt   - the gene set whose addition yields the largest control-energy
               reduction (legacy code: IER_CS; Table S1 abbrev: MER_CS).
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Dict, List, Set, Tuple

import networkx as nx
import numpy as np
from sympy import Matrix

from .io_utils import get_network_with_edge_dict, read_gene_list, write_gene_list


# --------------------------------------------------------------------------- #
# helpers (verbatim from notebook)
# --------------------------------------------------------------------------- #
def _regulators_of(G: nx.DiGraph, D: Set[str]) -> Set[str]:
    """All upstream regulators of disease core (not in D itself)."""
    r: List[str] = []
    for node in D:
        if node not in G:
            continue
        for i in G.predecessors(node):
            r.append(i)
    return set(r).difference(D)


def _bigraph(G: nx.DiGraph, driver: Set[str], dis: Set[str]) -> nx.DiGraph:
    bG = nx.DiGraph()
    for i in driver:
        for j in dis:
            if j in G.successors(i):
                bG.add_edge(i, j)
                bG.add_node(i)
                bG.add_node(j)
    return bG


def _matrix(
    bG: nx.DiGraph,
    dis: Set[str],
    driver: Set[str],
    edge_dict: Dict[str, List[str]],
) -> np.ndarray:
    """Verbatim port of the optimised ``matrix`` from the notebook."""
    dic: dict = {}
    A = np.zeros(shape=[len(set(dis)), len(set(driver))])
    i = 0
    j = 0
    for node in bG.nodes():
        if node in driver:
            dic[node] = i
            i += 1
        if node in dis:
            dic[node] = j
            j += 1
    for g1 in set(driver):
        h = dic[g1]
        add = list(set(edge_dict.get(g1, [])).intersection(set(dis)))
        for item in add:
            A[dic[item]][h] = 1
    return A


def _max_lin_indep_group(A: np.ndarray) -> np.ndarray:
    """Maximum linearly independent group via sympy rref.

    Verbatim from the notebook.  This relies on sympy 1.4 / numpy 1.17
    behaviour: ``Matrix(A).rref()`` on a float64 matrix correctly zeroes
    the dependent rows.  Use the conda env pinned in ``environment.yml``
    (or ``requirements.txt``) to reproduce.

    NOTE: newer sympy (>= ~1.5) leaves tiny float residues in dependent
    rows so ``(B == 0).all(1)`` misses them, the resulting A21 is rank
    deficient, and ``det(A21 @ A21.T)`` becomes 0 -> NaN energy.  The
    locked legacy env is the supported configuration.
    """
    M = Matrix(A).rref()[0]
    B = np.array(M)
    M1 = B[~(B == 0).all(1)]
    return np.mat(M1, "int")


def _control_energy(A21: np.ndarray) -> float:
    """E_u expectation, verbatim from the notebook."""
    M = A21.shape[0]
    N = A21.shape[1]
    asin = 0.0
    for i in range(M):
        A = np.delete(A21, i, 0)
        detQ = np.linalg.det(np.dot(A21, A21.T))
        detQkk = np.linalg.det(np.dot(A, A.T))
        asin += detQkk / detQ
    return 1.0 / 3.0 * N + M + 4 * asin


def _random_add(
    G: nx.DiGraph,
    bG: nx.DiGraph,
    regulators: Set[str],
    driver: Set[str],
    dis: Set[str],
    edge_dict: Dict[str, List[str]],
    rng: random.Random,
) -> Tuple[nx.DiGraph, List[str], List[str]]:
    """Verbatim port of ``random_add_node``."""
    candidate_drivers: List[str] = []
    for i in regulators:
        for j in dis:
            if (i, j) in G.edges():
                candidate_drivers.append(i)

    if len(dis) <= 100:
        nadd = round(len(dis) / 5)
    else:
        nadd = 20

    pool = sorted(set(candidate_drivers).difference(driver))
    if nadd > len(pool):
        nadd = len(pool)
    listadd = rng.sample(pool, nadd) if nadd > 0 else []

    drivernew = list(driver)
    for i in listadd:
        add = list(set(edge_dict.get(i, [])).intersection(set(dis)))
        if len(add) > 0:
            drivernew.append(i)
            bG.add_node(i)
            for item in add:
                bG.add_edge(i, item)
    return bG, list(set(drivernew)), listadd


# --------------------------------------------------------------------------- #
# top-level
# --------------------------------------------------------------------------- #
def run(
    network_path: str | Path,
    disease_core_path: str | Path,
    cp_ind_path: str | Path,
    cp_ext_path: str | Path,
    cp_drv_part1_path: str | Path,
    cp_drv_part2_path: str | Path,
    out_dir: str | Path,
    *,
    n_random: int = 100,
    seed: int | None = None,
) -> List[str]:
    """Compute CP-EFF.  Returns the gene list with smallest control energy."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    G, net_genes, edge_dict = get_network_with_edge_dict(network_path)
    disease_core = read_gene_list(disease_core_path)
    cp_ind = read_gene_list(cp_ind_path)
    cp_ext = read_gene_list(cp_ext_path)
    cp_drv_part1 = read_gene_list(cp_drv_part1_path)
    cp_drv_part2 = read_gene_list(cp_drv_part2_path)

    rng = random.Random(seed) if seed is not None else random.Random()

    # Sort set-derived collections so the same --seed gives reproducible
    # results across processes (Python hash randomisation otherwise makes
    # set iteration order vary, which propagates through _bigraph / _matrix
    # column ordering and ultimately changes the energy value).  Sorting
    # does not change which genes are selected, only the deterministic
    # tie-breaking inside numerical-rank operations.
    disnode = sorted(set(disease_core).intersection(net_genes))
    driver = sorted(set(cp_ind).union(cp_ext).union(cp_drv_part1).union(cp_drv_part2))
    regulator_nodes = sorted(_regulators_of(G, set(disnode)))

    bG = _bigraph(G, driver, disnode)
    A1 = _matrix(bG, disnode, driver, edge_dict)
    if np.linalg.matrix_rank(A1) == len(set(disnode)):
        e_orig = _control_energy(A1)
    else:
        A21 = _max_lin_indep_group(A1)
        e_orig = _control_energy(A21)

    rand = n_random
    result = np.zeros(rand + 1)
    result[0] = e_orig
    min_e = 0.0
    min_gene: List[str] = []

    for j in range(rand):
        bGnew, drivernew, listadd = _random_add(
            G, bG, regulator_nodes, driver, disnode, edge_dict, rng
        )
        A1new = _matrix(bGnew, disnode, set(drivernew), edge_dict)
        if np.linalg.matrix_rank(A1) == len(set(disnode)):
            b = _control_energy(A1new)
        else:
            A21new = _max_lin_indep_group(A1new)
            b = _control_energy(A21new)
        result[j + 1] = b
        if b < result[j]:
            min_gene_local: List[str] = []
            min_e = b
            if min_e < e_orig:
                min_gene_local = listadd
            min_gene = min_gene_local

    write_gene_list(min_gene, out_dir / "CP-EFF.txt")
    print(
        f"[Step 8] CP-EFF: original energy {e_orig:.6f}, "
        f"minimal energy {min_e:.6f}, {len(min_gene)} genes"
    )
    return min_gene
