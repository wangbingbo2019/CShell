"""Step 6 - obtain CP-DRV (part 1) and CP-SEN.

1:1 port of `[3] Obtain DCS_CS_1 or SOS_CS.ipynb`.

Reads the ``.ControlCentrality`` files produced by Yang-Yu Liu's external
CalControlCentrality tool (Step 5, see THIRD_PARTY.md), keeps only nodes NOT
already in the disease core, sorts by control-centrality value descending,
and selects the top 1% as the CP-DRV (part 1) / CP-SEN gene lists.

Output files (paper terminology)
--------------------------------
CP-DRV_part1.txt   - top 1% of nodes by control-subspace dimension
                     (legacy code: DCS_CS_1; Table S1 abbrev: CLS_CS part 1)
CP-SEN.txt         - top 1% of nodes by observation-subspace dimension
                     (legacy code: SOS_CS; Table S1 abbrev: OLS_CS)

The pipeline aggregates ``CP-DRV_part1.txt`` (this step) and
``CP-DRV_part2.txt`` (Step 7) into the final ``CP-DRV.txt``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

from .io_utils import read_gene_list, write_gene_list


def _read_sorted_subspace(
    cc_file: str | Path,
    disease_core: List[str],
) -> Dict[str, int]:
    """Return ``{gene: subspace_size}`` sorted descending, disease genes excluded.

    Verbatim port of ``read_input`` from the notebook.
    """
    dis = set(disease_core)

    # first pass: collect node names
    g: List[str] = []
    with open(cc_file, "r") as data:
        next(data)  # skip header
        for line in data:
            pas = line.strip().split()
            if not pas:
                continue
            g.append(pas[0])

    g1 = set(g).difference(dis)

    # second pass: collect values for non-disease genes
    genesub: Dict[str, int] = {}
    with open(cc_file, "r") as data:
        next(data)
        for line in data:
            pas = line.strip().split()
            if not pas:
                continue
            gene = pas[0]
            value = int(pas[1])
            if gene in g1:
                genesub[gene] = value

    sorted_items = sorted(genesub.items(), key=lambda kv: kv[1], reverse=True)
    return dict(sorted_items)


def _top_one_percent(sorted_dict: Dict[str, int]) -> List[str]:
    """Verbatim port of ``rank_choose``: keep ties at the 1% cutoff."""
    n = round(len(sorted_dict) * 0.01)
    values = list(sorted_dict.values())
    keys = list(sorted_dict.keys())
    if not values:
        return []
    a = values[n]
    return [keys[i] for i in range(len(sorted_dict)) if values[i] >= a]


def run(
    regulation_cc_file: str | Path,
    observation_cc_file: str | Path,
    disease_core_path: str | Path,
    out_dir: str | Path,
) -> Tuple[List[str], List[str]]:
    """Compute CP-DRV (part 1) and CP-SEN.

    Returns
    -------
    (cp_drv_part1_genes, cp_sen_genes)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    disease_core = read_gene_list(disease_core_path)

    sorted_subspace = _read_sorted_subspace(regulation_cc_file, disease_core)
    cp_drv_part1 = _top_one_percent(sorted_subspace)

    sorted_observ = _read_sorted_subspace(observation_cc_file, disease_core)
    cp_sen = _top_one_percent(sorted_observ)

    write_gene_list(cp_drv_part1, out_dir / "CP-DRV_part1.txt")
    write_gene_list(cp_sen, out_dir / "CP-SEN.txt")

    print(f"[Step 6] CP-DRV part1: {len(cp_drv_part1)} genes")
    print(f"[Step 6] CP-SEN:       {len(cp_sen)} genes")
    return cp_drv_part1, cp_sen
