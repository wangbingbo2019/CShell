# CShell pipeline overview

This document maps the eight steps in the paper to the code in this
repository, identifies which steps are *third-party* and which are
*ours*, and details the algorithm of each step that we implement.

For terminology, please refer back to the README.  In short:

* The paper's upper labels are **CP-IND**, **CP-EXT**, **CP-DRV**,
  **CP-SEN**, **CP-EFF** (used everywhere in this repo).
* Earlier code revisions and the supplementary Table S1 of the paper
  use the abbreviations ISR_CS, ESR_CS, CLS_CS (= DCS_CS_1 + DCS_CS_2),
  OLS_CS (= SOS_CS) and MER_CS (= IER_CS).  These are kept only as
  cross-references in `cshell/__init__.py`.

---

## Step 1 - Global controllability analysis (third-party)

  * Tool: **ControllabilityAnalysis** by Yang-Yu Liu (Nature 2011).
  * Status: NOT in this repository.  See [THIRD_PARTY.md](../THIRD_PARTY.md).
  * Inputs: a directed network edgelist.
  * Outputs: `*.nodetype`, `*.linktype`, `*.output`, `*.nodemap`.

For the default network we ship under `data/network/`, the outputs
`Directed_network.txt.nodetype` and `Directed_network.txt.linktype`
have already been generated and committed - so you do **not** need
to run Step 1 to reproduce the asthma example.

---

## Step 2 - Global Indispensable / Extraordinary node sets (precomputed)

This step was implemented by `[0] Control signature.ipynb` in earlier
revisions.  For the default network we ship the resulting
`Indispensable.txt` and `Extraordinary.txt` under `data/network/`,
so this step does NOT need to be re-run for the asthma example.

If you change the network, regenerate these files by re-running Step 1
and then the original control-signature logic.

---

## Step 3 - CP-IND and CP-EXT (ours)

  * Module: `cshell/step3_cp_ind_ext.py`
  * Source: 1:1 port of `[1] Obtain ISR_CS and ESR_CS.ipynb`.

Algorithm:

  1. Take the first-order neighbours `N` of the disease core `D`
     (both successors and predecessors of every gene in `D`,
     minus `D` itself).
  2. Intersect `N` with `Indispensable.txt` -> candidate set `I_n`.
     Intersect `N` with `Extraordinary.txt` -> candidate set `E_n`.
  3. For every candidate gene `g` in `I_n`, compute the
     **regulation ability** as `|successors(g) intersect D|`.
  4. Build a 100-sample random null model: for each iteration `i`
     sample a random subset `R_i` of the network of size `|D|` and
     recompute the regulation ability against `R_i`.
  5. The **z-score** of `g` is `(observed - mean(null)) / std(null)`.
  6. Keep only genes with z-score >= 2.33 -> `CP-IND` (legacy: `ISR_CS`).
  7. Repeat 3-6 with the candidate set `E_n` -> `CP-EXT` (legacy: `ESR_CS`).

The threshold 2.33 is the original notebook value, equivalent to a
one-tailed p-value of approximately 0.01.

---

## Step 4 - Regulation / observation subgraphs (ours)

  * Module: `cshell/step4_subgraphs.py`
  * Source: 1:1 port of `[2] Construct regulation subgraph and observation subgraph.ipynb`.

Builds two directed subgraphs that act as the input to the third-party
Step 5:

  * `regulation_subgraph.txt` - all edges from `(D u N)` *into* `D`.
  * `observation_subgraph.txt` - all edges from `D` *into* `(D u N)`.

The output filenames are kept as `regulation_subgraph.txt` and
`observation_subgraph.txt` so they remain drop-in compatible with
CalControlCentrality (which produces a sibling `*.ControlCentrality`
file).

---

## Step 5 - Control centrality (third-party)

  * Tool: **CalControlCentrality v1.1** by Yang-Yu Liu (PLoS ONE 2012).
  * Status: NOT in this repository.  See [THIRD_PARTY.md](../THIRD_PARTY.md).
  * Inputs: `regulation_subgraph.txt`, `observation_subgraph.txt`.
  * Outputs: `*.ControlCentrality` files containing per-node controllable
    / observable subspace dimensions.

For the asthma example we ship the resulting `*.ControlCentrality`
files under `data/examples/asthma/`, so the example reproduces without
having to install GLPK and compile the C++ tool.

---

## Step 6 - CP-DRV (part 1) and CP-SEN (ours)

  * Module: `cshell/step6_cp_drv_sen.py`
  * Source: 1:1 port of `[3] Obtain DCS_CS_1 or SOS_CS.ipynb`.

Algorithm:

  1. Read the controllable-subspace `.ControlCentrality` file produced
     by Step 5 on the regulation subgraph.
  2. Drop genes that are already in the disease core.
  3. Sort the remaining genes by controllable-subspace dimension
     (descending) and keep the top 1% -> `CP-DRV` part 1
     (legacy: `DCS_CS_1`; Table S1: `CLS_CS` part 1).
  4. Repeat 1-3 with the observability `.ControlCentrality` file ->
     `CP-SEN` (legacy: `SOS_CS`; Table S1: `OLS_CS`).

The "top 1%" is implemented as the notebook does: rank the values, find
the 1%-th value, and keep all genes whose value is `>=` that cutoff
(so any ties at the boundary are kept together).

---

## Step 7 - CP-DRV (part 2) (ours)

  * Module: `cshell/step7_cp_drv_extra.py`
  * Source: 1:1 port of `[4] obtain DCS_CS_2(Increase control range).ipynb`.

Algorithm:

  1. Seed `driver = CP-IND u CP-EXT` (the union of the Step 3 outputs).
  2. Build the bipartite incidence matrix `B[disease_core, driver]`
     where `B[d, g] = 1` iff `(g, d)` is an edge of the global network.
  3. Compute `r0 = rank(B)`.
  4. Sort all *non-driver, non-disease-core* nodes by their
     controllable-subspace dimension (descending).
  5. For each candidate gene `c`, compute `rank(B u {c})`.  If it is
     strictly greater than `r0`, add `c` to the driver set, update
     `r0`, and continue.  Stop when `r0 == |disease_core|` (full
     structural controllability) or no candidate increases the rank.
  6. The genes that were added in step 5 form `CP-DRV` part 2
     (legacy: `DCS_CS_2`; Table S1: `CLS_CS` part 2).

The pipeline aggregates `CP-DRV_part1.txt` and `CP-DRV_part2.txt` into
the final `CP-DRV.txt` (the union, sorted alphabetically).

---

## Step 8 - CP-EFF (ours)

  * Module: `cshell/step8_cp_eff.py`
  * Source: 1:1 port of `[5] Obtain IER_CS.ipynb`.

Algorithm (random search):

  1. Seed `driver = CP-IND u CP-EXT u CP-DRV (part1+part2)`.
  2. Build the bipartite incidence matrix `B` of `driver` vs.
     `disease_core`; if it is rank-deficient, replace it by its maximum
     linearly independent group of rows.
  3. Compute the **expected control energy** `E_u = N/3 + M + 4 * S`,
     where `S = sum_{i=1..M} det(Q_kk) / det(Q)`,
     `Q = B B^T`, and `Q_kk` is `Q` with row `i` (and the corresponding
     row of `B`) removed.  This is the closed-form expression used in
     the notebook; see the paper supplementary for the derivation.
  4. For 100 random iterations:
       a. Randomly sample a small subset of *upstream regulators* of
          the disease core (5 if `|D| <= 100` else 20 nodes - the
          notebook formula `nadd = round(|D|/5)` capped at 20).
       b. Recompute `E_u` after adding the sample to `driver`.
       c. Keep the **smallest** `E_u` seen so far and remember the
          random sample that produced it.
  5. The recorded random sample is `CP-EFF` (legacy: `IER_CS`;
     Table S1: `MER_CS`).

The default 100-iteration random search matches the notebook.  Since
this step is non-deterministic, the gene list will differ between runs
unless you pass `--seed` to `scripts/run_pipeline.py`.

---

## End-to-end orchestration

  * Module: `cshell/pipeline.py`
  * CLI: `scripts/run_pipeline.py`

The pipeline calls Step 3 -> Step 4 -> Step 6 -> Step 7 -> Step 8 in
order, then writes the unified `cshell_genes.tsv` with two columns:

```
gene<TAB>cp_class
```

`cp_class` takes one of `CP-IND`, `CP-EXT`, `CP-DRV`, `CP-SEN`,
`CP-EFF`.  A gene that satisfies multiple CP definitions appears on
multiple rows (one per class).  This format is the canonical answer
to the *"which CP class did this gene come from?"* question and was
chosen for direct compatibility with the paper's tables.
