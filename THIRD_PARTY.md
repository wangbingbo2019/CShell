# Third-party tools and attribution

The CShell pipeline depends on **two C++ packages developed by Yang-Yu
Liu and colleagues** at the Channing Division of Network Medicine,
Harvard Medical School.  These tools implement the underlying network
controllability and control-centrality theories and are **not part of
this repository**.  They are not redistributed here because their
license explicitly forbids redistribution.

This document explains:

  1. Which steps of the paper they cover.
  2. How to obtain them.
  3. Their license restrictions.
  4. How their outputs flow into CShell's Python pipeline.

================================================================================

## 1. The two third-party packages

| Package | Paper step | Used to compute | Reference |
|---------|------------|-----------------|-----------|
| ControllabilityAnalysis | Step 1 | Per-node and per-edge controllability classes (`*.nodetype`, `*.linktype`). These feed paper Step 2 (the `[0] Control signature.ipynb` notebook in older revisions of this code) which enumerates the global *Indispensable* and *Extraordinary* node sets. | Liu Y-Y, Slotine J-J, Barabasi A-L. *Controllability of complex networks*. **Nature** 473, 167-173 (12 May 2011). |
| CalControlCentrality v1.1 | Step 5 | The control centrality (controllable subspace dimension and observable subspace dimension) for every node of the regulation / observation subgraph built around the disease core. Output: `*.ControlCentrality`. | Liu Y-Y, Slotine J-J, Barabasi A-L. *Control Centrality and Hierarchical Structure in Complex Networks*. **PLoS ONE** 7(9): e44459 (2012). doi:10.1371/journal.pone.0044459 |

Both tools were obtained by us as a *distribution package* directly from
Yang-Yu Liu (the original author).  This is consistent with the original
software-usage note that documented our setup: *"We get this distribution
package from Yang-Yu Liu."*  We do **not** claim authorship over them,
nor do we re-host their source.

================================================================================

## 2. License restrictions (verbatim from the original READMEs)

> The permission to use this package is granted by the author Yang-Yu
> Liu, only on the following 5 conditions:
>
> 1. The package remains at your institution and should not be
>    published, distributed, or otherwise transferred or made available
>    to people other than institution employees and students involved in
>    research under your supervision.
>
> 2. The package should be used by you and/or your institution solely
>    for non-commercial purposes, except with express permission from
>    the author.
>
> 3. You may provide the author with feedback on the use of the package
>    in your research, and that the author is permitted to use any
>    information you provide in making changes to the package.
>
> 4. Any risk associated with using the package at your institution is
>    with you and your institution.
>
> 5. The package should be cited in any publication(s) reporting on
>    data obtained by using it.

Because of clause 1, this GitHub repository **must not** include the
source code, compiled binaries, or compressed archives of either
package.  The `.gitignore` file is configured to keep them out by
default.

================================================================================

## 3. How to obtain the third-party packages

Contact the corresponding author (Yang-Yu Liu) directly:

    Address: Channing Division of Network Medicine,
             181 Longwood Avenue, Boston, MA 02115
    Email:   yyl@channing.harvard.edu

================================================================================

## 4. Data flow (how the third-party outputs feed CShell)

```
+-------------------------+         +-------------------------------+
|  Directed_network.txt   |  ---->  |  Step 1: ControllabilityAnalysis |
+-------------------------+         |  (third-party C++)            |
                                    +---------------+---------------+
                                                    |
                                                    v
                                  +---------------------------------+
                                  |  Directed_network.txt.nodetype   |
                                  |  Directed_network.txt.linktype   |
                                  +---------------+---------------+
                                                    |
                                                    v
                                  +---------------------------------+
                                  |  Step 2: Control signature       |
                                  |  ([0] Control signature.ipynb)   |
                                  |  -- pre-computed for the default |
                                  |     network and shipped under    |
                                  |     data/network/                |
                                  +---------------+---------------+
                                                    |
                                                    v
                                  +---------------------------------+
                                  |  Indispensable.txt               |
                                  |  Extraordinary.txt               |
                                  +---------------+---------------+
                                                    |
                                                    v   (network-level)
+-------------------------+         +---------------------------------+
|  disease_core.txt       |  -+-->  |  CShell Step 3 (CP-IND, CP-EXT)  |
+-------------------------+   |     +---------------+---------------+
                              |                     |
                              |                     v
                              +-->  +---------------------------------+
                                    |  CShell Step 4 (subgraphs)       |
                                    +---------------+---------------+
                                                    |
                                                    v
                                  +---------------------------------+
                                  |  regulation_subgraph.txt         |
                                  |  observation_subgraph.txt        |
                                  +---------------+---------------+
                                                    |
                                                    v
                                  +---------------------------------+
                                  |  Step 5: CalControlCentrality    |
                                  |  (third-party C++)               |
                                  +---------------+---------------+
                                                    |
                                                    v
                                  +---------------------------------+
                                  |  *.ControlCentrality             |
                                  +---------------+---------------+
                                                    |
                                                    v
                                  +---------------------------------+
                                  |  CShell Step 6 (CP-DRV part1,    |
                                  |                 CP-SEN)          |
                                  |  CShell Step 7 (CP-DRV part2)    |
                                  |  CShell Step 8 (CP-EFF)          |
                                  +---------------+---------------+
                                                    |
                                                    v
                                  +---------------------------------+
                                  |  cshell_genes.tsv  (unified)     |
                                  |  CP-IND.txt, CP-EXT.txt, ...     |
                                  +---------------------------------+
```

================================================================================

## 5. Re-running Step 1 / Step 5 yourself

If you change the network you must re-run Step 1+2 once for the new
network (this regenerates `*.nodetype`, `*.linktype`, `Indispensable.txt`
and `Extraordinary.txt`).  If you change the disease core you must
re-run Step 5 once for the new disease (this regenerates the two
`*.ControlCentrality` files).

The third-party READMEs explain installation and usage:

    Step 1 (ControllabilityAnalysis):
        ./Parse                ./Directed_network.txt
        ./ControllabilityAnalysis ./Directed_network.txt
        # produces .nodetype, .linktype, .nodemap, .output

    Step 5 (CalControlCentrality, requires GLPK):
        ./CalControlCentrality . regulation_subgraph.txt
        ./CalControlCentrality . observation_subgraph.txt
        # produces ./data/regulation_subgraph.txt.ControlCentrality
        # and     ./data/observation_subgraph.txt.ControlCentrality

Then place the resulting files at the paths CShell expects (see
`README.md`, "Inputs") and re-run `scripts/run_pipeline.py`.

The two third-party READMEs are reproduced **only inside the original
distribution packages**; we do not include them here.
