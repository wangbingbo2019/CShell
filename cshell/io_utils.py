"""I/O helpers shared by all pipeline steps.

Functions are direct extracts from the original notebooks (`get_network`,
`readinput`, `read_input`); they are unified here so each step does not
re-define them.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import networkx as nx
import pandas as pd


# --------------------------------------------------------------------------- #
# network
# --------------------------------------------------------------------------- #
def get_network(network_path: str | Path) -> Tuple[nx.DiGraph, List[str]]:
    """Read a 2-column whitespace-separated edgelist into a directed graph.

    Returns the graph and its node list (in insertion order).
    """
    network_path = Path(network_path)
    G = nx.DiGraph()
    with open(network_path, "r") as f:
        for line in f:
            pas = line.strip().split()
            if len(pas) < 2:
                continue
            G.add_node(pas[0])
            G.add_node(pas[1])
            G.add_edge(pas[0], pas[1])
    return G, list(G.nodes())


def get_network_with_edge_dict(
    network_path: str | Path,
) -> Tuple[nx.DiGraph, List[str], Dict[str, List[str]]]:
    """Same as :func:`get_network` but also returns an out-neighbour dict.

    Used by the IER (CP-EFF) step which needs O(1) successor lookup.
    """
    network_path = Path(network_path)
    G, node_list = get_network(network_path)
    all_edge = pd.read_table(network_path, header=None)
    edge_dict: Dict[str, List[str]] = {}
    for src in set(all_edge[0]):
        edge_dict[src] = []
    for i in range(len(all_edge)):
        edge_dict[all_edge.iloc[i, 0]].append(all_edge.iloc[i, 1])
    return G, node_list, edge_dict


# --------------------------------------------------------------------------- #
# gene lists
# --------------------------------------------------------------------------- #
def read_gene_list(path: str | Path) -> List[str]:
    """Read a one-gene-per-line text file into a list (preserves order)."""
    genes: List[str] = []
    with open(path, "r") as f:
        for line in f:
            n = line.strip().split()
            if not n:
                continue
            genes.append(n[0])
    return genes


def write_gene_list(genes: Iterable[str], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for g in genes:
            f.write(f"{g}\n")
