"""Step 7 - obtain CP-DRV (part 2) - "Increase control range".

1:1 port of `[4] obtain DCS_CS_2(Increase control range).ipynb`.

Greedy algorithm: starting from driver = CP-IND u CP-EXT, scan candidate
genes (sorted by controllable subspace) and add any that strictly increase
the rank of the driver-vs-disease bipartite incidence matrix; stop when
the rank reaches |disease_core| (full controllability) or no more rank
improvement is possible.

Output file (paper terminology)
-------------------------------
CP-DRV_part2.txt   - bipartite-rank fillers added on top of CP-IND u CP-EXT
                     (legacy code: DCS_CS_2; Table S1 abbrev: CLS_CS part 2)

The pipeline aggregates ``CP-DRV_part1.txt`` (Step 6) and ``CP-DRV_part2.txt``
(this step) into the final ``CP-DRV.txt``.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Set, Tuple

import networkx as nx
import numpy as np

from .io_utils import get_network, read_gene_list, write_gene_list


def _sort_subspace_excluding(
    cc_file: str | Path,
    disease_list: List[str],
    driver_list: List[str],
) -> List[str]:
    """Verbatim port of ``sort_subspace`` from the notebook."""
    gene_in_sub: dict = {}
    with open(cc_file, "r") as data:
        next(data)  # header
        for line in data:
            pas = line.strip().split()
            if not pas:
                continue
            if pas[0] not in disease_list and pas[0] not in driver_list:
                gene_in_sub[pas[0]] = int(pas[1])

    sorted_items = sorted(gene_in_sub.items(), key=lambda kv: kv[1], reverse=True)
    return list(dict(sorted_items).keys())


def _bigraph(
    G: nx.DiGraph, driver: List[str], dis: Set[str]
) -> Tuple[nx.DiGraph, Set[str], Set[str]]:
    """Verbatim port of ``bigraph`` from the notebook."""
    bG = nx.DiGraph()
    dis_in_bg: List[str] = []
    driver_in_bg: List[str] = []
    for i in driver:
        for j in dis:
            if j in G.successors(i):
                bG.add_edge(i, j)
                bG.add_node(i)
                bG.add_node(j)
                driver_in_bg.append(i)
                dis_in_bg.append(j)
    return bG, set(dis_in_bg), set(driver_in_bg)


def _matrix(bG: nx.DiGraph, driver_in_g: Set[str], dis_in_g: Set[str]) -> np.ndarray:
    """Verbatim port of ``matrix`` from the notebook."""
    dic: dict = {}
    A = np.zeros(shape=[len(set(dis_in_g)), len(set(driver_in_g))])

    i = 0
    j = 0
    for node in bG.nodes():
        if node in driver_in_g:
            dic[node] = i
            i += 1
    for node in bG.nodes():
        if node in dis_in_g:
            dic[node] = j
            j += 1
    for g1 in driver_in_g:
        h = dic[g1]
        for g2 in dis_in_g:
            k = dic[g2]
            if (g1, g2) in bG.edges():
                A[k][h] = 1
    return A


def _full_rank_grow(
    G: nx.DiGraph,
    driver: List[str],
    dis: Set[str],
    candidates: List[str],
) -> Tuple[int, List[str]]:
    """Verbatim port of ``fullrank`` from the notebook (greedy rank-growing)."""
    originalbG, dis_in_g, driver_in_g = _bigraph(G, driver, dis)
    originalA = _matrix(originalbG, driver_in_g, dis_in_g)
    originalr = np.linalg.matrix_rank(originalA)

    geneadd: List[str] = []
    newr = originalr
    for gene in candidates:
        if originalr == len(dis):
            newr = originalr
            break
        driver1 = set(driver).union({gene})
        newbG, dis_in_bg, driver_in_bg = _bigraph(G, list(driver1), dis)
        newA = _matrix(newbG, driver_in_bg, dis_in_bg)
        newr = np.linalg.matrix_rank(newA)

        if newr == len(dis) or len(geneadd) == len(dis):
            driver.append(gene)
            geneadd.append(gene)
            break
        elif newr > originalr:
            driver.append(gene)
            geneadd.append(gene)
            newbG1, dis_in_g, driver_in_g = _bigraph(G, driver, dis)
            newA1 = _matrix(newbG1, driver_in_g, dis_in_bg)
            originalr = np.linalg.matrix_rank(newA1)
        else:
            continue

    print(
        f"[Step 7] original rank {np.linalg.matrix_rank(originalA)} "
        f"full rank {len(dis)} new rank {newr} added {len(geneadd)} {geneadd}"
    )
    return newr, geneadd


def run(
    network_path: str | Path,
    disease_core_path: str | Path,
    cp_ind_path: str | Path,
    cp_ext_path: str | Path,
    regulation_cc_file: str | Path,
    out_dir: str | Path,
) -> List[str]:
    """Compute CP-DRV (part 2). Returns added gene list.

    The driver set seed is ``CP-IND u CP-EXT``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    G, net_genes = get_network(network_path)
    disease_core = read_gene_list(disease_core_path)
    cp_ind = read_gene_list(cp_ind_path)
    cp_ext = read_gene_list(cp_ext_path)

    disease = set(net_genes).intersection(disease_core)
    driver = list(set(cp_ind).union(cp_ext))
    driver.sort()

    candidates = _sort_subspace_excluding(regulation_cc_file, list(disease), driver)
    _, geneadd = _full_rank_grow(G, driver, disease, candidates)

    write_gene_list(geneadd, out_dir / "CP-DRV_part2.txt")
    print(f"[Step 7] CP-DRV part2: {len(geneadd)} genes")
    return geneadd
