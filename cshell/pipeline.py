"""End-to-end CShell pipeline.

Given a disease core gene list, a directed network, the pre-computed global
indispensable / extraordinary node lists, and the two ``.ControlCentrality``
files (Step 5 of the paper - run externally with Yang-Yu Liu's
CalControlCentrality, see THIRD_PARTY.md), this orchestrates Step 3, 4, 6, 7,
8 and writes a single TSV ``cshell_genes.tsv`` with two columns:

    gene<TAB>cp_class

Where cp_class is one of CP-IND, CP-EXT, CP-DRV, CP-SEN, CP-EFF.  A gene
selected by more than one class is emitted on multiple rows, one per class.

Per-class intermediate files (paper terminology, also written to ``out_dir``):

    CP-IND.txt
    CP-EXT.txt
    CP-DRV_part1.txt + CP-DRV_part2.txt   --merged into-->   CP-DRV.txt
    CP-SEN.txt
    CP-EFF.txt
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from . import (
    step3_cp_ind_ext,
    step4_subgraphs,
    step6_cp_drv_sen,
    step7_cp_drv_extra,
    step8_cp_eff,
)
from .io_utils import write_gene_list


def _write_cshell_tsv(
    classes: Dict[str, Iterable[str]],
    out_path: Path,
) -> int:
    """Write the unified ``gene\tcp_class`` long-format TSV.

    Order: by CP class (paper order), then by gene name within a class.
    Returns total row count.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows: List[Tuple[str, str]] = []
    for cp in ("CP-IND", "CP-EXT", "CP-DRV", "CP-SEN", "CP-EFF"):
        for gene in sorted(set(classes.get(cp, []))):
            rows.append((gene, cp))
    with open(out_path, "w") as f:
        f.write("gene\tcp_class\n")
        for gene, cp in rows:
            f.write(f"{gene}\t{cp}\n")
    return len(rows)


def run_full_pipeline(
    network_path: str | Path,
    disease_core_path: str | Path,
    indispensable_path: str | Path,
    extraordinary_path: str | Path,
    regulation_cc_file: str | Path,
    observation_cc_file: str | Path,
    out_dir: str | Path,
    *,
    z_threshold: float = 2.33,
    n_random_step3: int = 100,
    n_random_step8: int = 100,
    seed: int | None = None,
) -> Path:
    """Run all CShell steps and write ``cshell_genes.tsv`` under ``out_dir``.

    Returns the path to the unified TSV.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Step 3: CP-IND, CP-EXT
    cp_ind, cp_ext = step3_cp_ind_ext.run(
        network_path,
        disease_core_path,
        indispensable_path,
        extraordinary_path,
        out_dir,
        z_threshold=z_threshold,
        n_random=n_random_step3,
        seed=seed,
    )

    # Step 4: regulation/observation subgraphs (kept on disk so the output
    # directory remains self-contained for inspection or for re-running
    # the third-party CalControlCentrality).
    step4_subgraphs.run(network_path, disease_core_path, out_dir)

    # Step 6: CP-DRV (part 1), CP-SEN
    cp_drv_part1, cp_sen = step6_cp_drv_sen.run(
        regulation_cc_file,
        observation_cc_file,
        disease_core_path,
        out_dir,
    )

    # Step 7: CP-DRV (part 2)
    cp_drv_part2 = step7_cp_drv_extra.run(
        network_path,
        disease_core_path,
        out_dir / "CP-IND.txt",
        out_dir / "CP-EXT.txt",
        regulation_cc_file,
        out_dir,
    )

    # CP-DRV = union of part 1 and part 2
    cp_drv = sorted(set(cp_drv_part1).union(cp_drv_part2))
    write_gene_list(cp_drv, out_dir / "CP-DRV.txt")

    # Step 8: CP-EFF
    cp_eff = step8_cp_eff.run(
        network_path,
        disease_core_path,
        out_dir / "CP-IND.txt",
        out_dir / "CP-EXT.txt",
        out_dir / "CP-DRV_part1.txt",
        out_dir / "CP-DRV_part2.txt",
        out_dir,
        n_random=n_random_step8,
        seed=seed,
    )

    classes: Dict[str, Iterable[str]] = {
        "CP-IND": cp_ind,
        "CP-EXT": cp_ext,
        "CP-DRV": cp_drv,
        "CP-SEN": cp_sen,
        "CP-EFF": cp_eff,
    }
    tsv_path = out_dir / "cshell_genes.tsv"
    n_rows = _write_cshell_tsv(classes, tsv_path)

    print()
    print(f"[Pipeline] CShell TSV written to {tsv_path} ({n_rows} rows)")
    return tsv_path
