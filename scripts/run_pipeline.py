"""CLI entry point: run the full CShell pipeline from the command line.

Example
-------
    python scripts/run_pipeline.py \\
        --network               data/network/Directed_network.txt \\
        --disease-core          data/examples/asthma/disease_core.txt \\
        --indispensable         data/network/Indispensable.txt \\
        --extraordinary         data/network/Extraordinary.txt \\
        --regulation-cc         data/examples/asthma/regulation_subgraph.txt.ControlCentrality \\
        --observation-cc        data/examples/asthma/observation_subgraph.txt.ControlCentrality \\
        --out-dir               results/asthma

The two ``--*-cc`` files must be produced beforehand by the third-party
CalControlCentrality tool (see THIRD_PARTY.md).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# allow running this script as ``python scripts/run_pipeline.py`` without install
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cshell.pipeline import run_full_pipeline  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="Run the full CShell pipeline (Steps 3, 4, 6, 7, 8) and "
                    "label every output gene as CP-IND/CP-EXT/CP-DRV/CP-SEN/CP-EFF.",
    )
    p.add_argument("--network", required=True, type=Path,
                   help="Directed network edgelist (2 columns).")
    p.add_argument("--disease-core", required=True, type=Path,
                   help="Disease core gene list (one gene per line).")
    p.add_argument("--indispensable", required=True, type=Path,
                   help="Pre-computed Indispensable.txt for the network "
                        "(produced by paper Step 2 from .nodetype/.linktype).")
    p.add_argument("--extraordinary", required=True, type=Path,
                   help="Pre-computed Extraordinary.txt for the network.")
    p.add_argument("--regulation-cc", required=True, type=Path,
                   help="*.ControlCentrality file for the regulation subgraph "
                        "(produced by the third-party CalControlCentrality, "
                        "Step 5 - see THIRD_PARTY.md).")
    p.add_argument("--observation-cc", required=True, type=Path,
                   help="*.ControlCentrality file for the observation subgraph.")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Output directory.")
    p.add_argument("--z-threshold", type=float, default=2.33,
                   help="Regulation-ability z-score threshold for CP-IND / CP-EXT "
                        "(default: 2.33, as in the paper).")
    p.add_argument("--n-random-step3", type=int, default=100,
                   help="Number of random null-model samples for Step 3 "
                        "(CP-IND, CP-EXT) (default: 100).")
    p.add_argument("--n-random-step8", type=int, default=100,
                   help="Number of random search iterations for Step 8 "
                        "(CP-EFF) (default: 100).")
    p.add_argument("--seed", type=int, default=None,
                   help="Optional random seed for reproducibility (default: not "
                        "seeded, matching the original notebooks).")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    run_full_pipeline(
        network_path=args.network,
        disease_core_path=args.disease_core,
        indispensable_path=args.indispensable,
        extraordinary_path=args.extraordinary,
        regulation_cc_file=args.regulation_cc,
        observation_cc_file=args.observation_cc,
        out_dir=args.out_dir,
        z_threshold=args.z_threshold,
        n_random_step3=args.n_random_step3,
        n_random_step8=args.n_random_step8,
        seed=args.seed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
