"""Step 4 - construct regulation and observation subgraphs.

1:1 port of `[2] Construct regulation subgraph and observation subgraph.ipynb`.

These subgraphs feed Yang-Yu Liu's external CalControlCentrality tool
(Step 5, NOT included here - see THIRD_PARTY.md).

The output files have the legacy names ``regulation_subgraph.txt`` and
``observation_subgraph.txt`` so they remain drop-in compatible with the
external tool's expected naming.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Set, Tuple

import networkx as nx

from .io_utils import get_network, read_gene_list
from .step3_cp_ind_ext import first_order_neighbors


def regulation_subgraph(G: nx.DiGraph, D: Set[str], N: Iterable[str]) -> nx.DiGraph:
    """Edges from (D u N) into D - "who regulates the disease core"."""
    rG = nx.DiGraph()
    for i in set(D).union(N):
        for j in D:
            if (i, j) in G.edges():
                rG.add_node(i)
                rG.add_node(j)
                rG.add_edge(i, j)
    return rG


def observation_subgraph(G: nx.DiGraph, D: Set[str], N: Iterable[str]) -> nx.DiGraph:
    """Edges from D into (D u N) - "what the disease core regulates"."""
    oG = nx.DiGraph()
    for i in set(D).union(N):
        for j in D:
            if (j, i) in G.edges():
                oG.add_node(i)
                oG.add_node(j)
                oG.add_edge(i, j)
    return oG


def _write_edgelist(graph: nx.DiGraph, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for (i, j) in graph.edges():
            f.write(f"{i}\t{j}\n")


def run(
    network_path: str | Path,
    disease_core_path: str | Path,
    out_dir: str | Path,
) -> Tuple[Path, Path]:
    """Write ``regulation_subgraph.txt`` and ``observation_subgraph.txt``.

    Returns the two output paths.  These files must be fed to
    ``CalControlCentrality`` (Step 5) by the user, see THIRD_PARTY.md.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    G, net_genes = get_network(network_path)
    disease_core = read_gene_list(disease_core_path)

    disease_module = set(net_genes).intersection(disease_core)
    neighbours = first_order_neighbors(G, disease_module)

    rG = regulation_subgraph(G, disease_module, neighbours)
    oG = observation_subgraph(G, disease_module, neighbours)

    reg_path = out_dir / "regulation_subgraph.txt"
    obs_path = out_dir / "observation_subgraph.txt"
    _write_edgelist(rG, reg_path)
    _write_edgelist(oG, obs_path)

    print(f"[Step 4] regulation_subgraph: {rG.number_of_edges()} edges -> {reg_path}")
    print(f"[Step 4] observation_subgraph: {oG.number_of_edges()} edges -> {obs_path}")
    return reg_path, obs_path
