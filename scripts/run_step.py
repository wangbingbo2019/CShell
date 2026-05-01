"""Run a single CShell step.

Useful for debugging or for partial reruns, e.g. when only the disease core
changed but the network's ``.ControlCentrality`` files are unchanged.

Example
-------
    python scripts/run_step.py step3 \\
        --network data/network/Directed_network.txt \\
        --disease-core data/examples/asthma/disease_core.txt \\
        --indispensable data/network/Indispensable.txt \\
        --extraordinary data/network/Extraordinary.txt \\
        --out-dir results/asthma
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cshell import (  # noqa: E402
    step3_cp_ind_ext,
    step4_subgraphs,
    step6_cp_drv_sen,
    step7_cp_drv_extra,
    step8_cp_eff,
)


STEPS = ("step3", "step4", "step6", "step7", "step8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="run_step.py")
    p.add_argument(
        "step", choices=STEPS,
        help="step3=CP-IND/CP-EXT, step4=subgraphs, "
             "step6=CP-DRV(part1)/CP-SEN, step7=CP-DRV(part2), step8=CP-EFF",
    )
    p.add_argument("--network", type=Path)
    p.add_argument("--disease-core", type=Path)
    p.add_argument("--indispensable", type=Path)
    p.add_argument("--extraordinary", type=Path)
    p.add_argument("--regulation-cc", type=Path)
    p.add_argument("--observation-cc", type=Path)
    p.add_argument("--cp-ind", type=Path, help="CP-IND.txt for steps 7 and 8.")
    p.add_argument("--cp-ext", type=Path, help="CP-EXT.txt for steps 7 and 8.")
    p.add_argument("--cp-drv-part1", type=Path, help="CP-DRV_part1.txt for step 8.")
    p.add_argument("--cp-drv-part2", type=Path, help="CP-DRV_part2.txt for step 8.")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--seed", type=int, default=None)
    return p.parse_args()


def _require(args: argparse.Namespace, *names: str) -> None:
    missing = [n for n in names if getattr(args, n.replace("-", "_")) is None]
    if missing:
        raise SystemExit(
            f"missing arguments for this step: {', '.join('--' + n for n in missing)}"
        )


def main() -> int:
    args = parse_args()
    if args.step == "step3":
        _require(args, "network", "disease-core", "indispensable", "extraordinary")
        step3_cp_ind_ext.run(
            args.network, args.disease_core, args.indispensable,
            args.extraordinary, args.out_dir, seed=args.seed,
        )
    elif args.step == "step4":
        _require(args, "network", "disease-core")
        step4_subgraphs.run(args.network, args.disease_core, args.out_dir)
    elif args.step == "step6":
        _require(args, "regulation-cc", "observation-cc", "disease-core")
        step6_cp_drv_sen.run(
            args.regulation_cc, args.observation_cc,
            args.disease_core, args.out_dir,
        )
    elif args.step == "step7":
        _require(args, "network", "disease-core", "cp-ind", "cp-ext", "regulation-cc")
        step7_cp_drv_extra.run(
            args.network, args.disease_core, args.cp_ind, args.cp_ext,
            args.regulation_cc, args.out_dir,
        )
    elif args.step == "step8":
        _require(
            args, "network", "disease-core",
            "cp-ind", "cp-ext", "cp-drv-part1", "cp-drv-part2",
        )
        step8_cp_eff.run(
            args.network, args.disease_core, args.cp_ind, args.cp_ext,
            args.cp_drv_part1, args.cp_drv_part2, args.out_dir, seed=args.seed,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
