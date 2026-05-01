"""Step 3 - obtain CP-IND and CP-EXT.

1:1 port of `[1] Obtain ISR_CS and ESR_CS.ipynb`.  The algorithm
(regulation-ability z-score against a 100-sample random null model,
threshold z >= 2.33) is preserved verbatim from the notebook - do NOT
modify the formula.

Output files (paper terminology)
--------------------------------
CP-IND.txt   - indispensable-derived regulators of disease core
               (legacy code variable: ISR_CS; Table S1 abbreviation: ISR_CS)
CP-EXT.txt   - extraordinary-derived regulators of disease core
               (legacy code variable: ESR_CS; Table S1 abbreviation: ESR_CS)
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import networkx as nx
import numpy as np

from .io_utils import get_network, read_gene_list, write_gene_list


# --------------------------------------------------------------------------- #
# helpers (verbatim from notebook)
# --------------------------------------------------------------------------- #
def first_order_neighbors(G: nx.DiGraph, D: Iterable[str]) -> Set[str]:
    """Both-direction first-order neighbours of ``D`` minus ``D`` itself."""
    D_set = set(D)
    neigh_1: List[str] = []
    for node in D_set:
        if node not in G:
            continue
        for i in G.successors(node):
            neigh_1.append(i)
        for i in G.predecessors(node):
            neigh_1.append(i)
    return set(neigh_1).difference(D_set)


def regulation_ability(
    G: nx.DiGraph,
    D: Set[str],
    candidates: List[str],
    net_genes: Sequence[str],
    n_random: int = 100,
    rng: random.Random | None = None,
) -> Dict[str, float]:
    """Compute regulation-ability z-score for every gene in ``candidates``.

    Verbatim from the notebook (`Regulation_Ability`).  ``n_random=100`` and
    the z-score normalisation `(observed - mean) / std` are exactly as in the
    original code.
    """
    if rng is None:
        rng = random.Random()  # not seeded -> matches notebook behaviour

    a = n_random
    Result = np.zeros((len(candidates)))
    Result1 = np.zeros((len(candidates), a + 1))

    # observed regulation strength: |successors(gene) intersect D|
    for gene in candidates:
        b = sum(1 for it in G.successors(gene) if it in D)
        Result[candidates.index(gene)] = b

    for i in range(len(candidates)):
        Result1[i, 0] = int(Result[i])

    # null model: 100 random gene sets the size of D
    D_size = len(D)
    for i in range(a):
        list1 = rng.sample(list(net_genes), D_size)
        random_set = set(list1)
        for gene in candidates:
            b = sum(1 for it in G.successors(gene) if it in random_set)
            Result1[candidates.index(gene), i + 1] = b

    meana = np.mean(Result1, axis=1)
    stda = np.std(Result1, axis=1)

    z_scores: Dict[str, float] = {}
    with np.errstate(invalid="ignore", divide="ignore"):
        for i in range(len(candidates)):
            # NaN propagates when stda==0; threshold filter naturally excludes
            z_scores[candidates[i]] = (Result[i] - meana[i]) / stda[i]
    return z_scores


def choose_above(scores: Dict[str, float], threshold: float) -> List[str]:
    """Keep genes with score >= ``threshold`` (NaN excluded)."""
    return [g for g, s in scores.items() if s >= threshold]


# --------------------------------------------------------------------------- #
# top-level entry
# --------------------------------------------------------------------------- #
def run(
    network_path: str | Path,
    disease_core_path: str | Path,
    indispensable_path: str | Path,
    extraordinary_path: str | Path,
    out_dir: str | Path,
    *,
    z_threshold: float = 2.33,
    n_random: int = 100,
    seed: int | None = None,
) -> Tuple[List[str], List[str]]:
    """Compute CP-IND and CP-EXT; write CP-IND.txt and CP-EXT.txt.

    Returns
    -------
    (cp_ind_genes, cp_ext_genes)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    G, net_genes = get_network(network_path)
    disease_core = read_gene_list(disease_core_path)
    indispensable = set(read_gene_list(indispensable_path))
    extraordinary = set(read_gene_list(extraordinary_path))

    rng = random.Random(seed) if seed is not None else None

    disease_module = set(net_genes).intersection(disease_core)
    neighbours = first_order_neighbors(G, disease_module)

    ind_in_neigh = list(neighbours.intersection(indispensable))
    ext_in_neigh = list(neighbours.intersection(extraordinary))

    z_ind = regulation_ability(G, disease_module, ind_in_neigh, net_genes, n_random, rng)
    z_ext = regulation_ability(G, disease_module, ext_in_neigh, net_genes, n_random, rng)

    cp_ind = choose_above(z_ind, z_threshold)
    cp_ext = choose_above(z_ext, z_threshold)

    write_gene_list(cp_ind, out_dir / "CP-IND.txt")
    write_gene_list(cp_ext, out_dir / "CP-EXT.txt")

    print(f"[Step 3] CP-IND: {len(cp_ind)} genes")
    print(f"[Step 3] CP-EXT: {len(cp_ext)} genes")
    return cp_ind, cp_ext
