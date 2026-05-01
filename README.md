# CShell - Control Shell of a Disease Core

This is the distribution package of CShell v1.0
2026/05/01

Authors: [PLEASE FILL: your name(s)]
Address: [PLEASE FILL: your institution]
Email:   [PLEASE FILL: your contact email]

CShell takes a disease *core* gene list and a directed gene-regulatory
network and returns the genes belonging to the five paper-defined
**control properties (CP)**:

| Label    | Concept (paper Table S1)               | Detection index           |
|----------|----------------------------------------|---------------------------|
| CP-IND   | Minimal driver set (indispensable)     | Z-score(i,d) > 2.33       |
| CP-EXT   | Control signature (extraordinary)      | Z-score(i,d) > 2.33       |
| CP-DRV   | Control subspace                       | top 1% (P-value < 0.01)   |
| CP-SEN   | Control observability                  | top 1% (P-value < 0.01)   |
| CP-EFF   | Control energy                         | Max(Δenergy)              |

Earlier code revisions used the abbreviations ISR_CS, ESR_CS,
DCS_CS_1/DCS_CS_2 (or CLS_CS), SOS_CS (or OLS_CS) and IER_CS (or MER_CS);
this repository uses **only** the upper paper labels CP-IND, CP-EXT,
CP-DRV, CP-SEN, CP-EFF in all output files and column names.


--------------------------------------------------------------------------
I. LICENSE
--------------------------------------------------------------------------
1. The Python source code in this repository (the directories `cshell/`
   and `scripts/` and the documentation we authored) is released under
   the MIT License - see [LICENSE](LICENSE).

2. CShell relies on two **third-party C++ packages** by Yang-Yu Liu
   (Channing Division of Network Medicine, Harvard Medical School):

      ControllabilityAnalysis      - paper Step 1
      CalControlCentrality v1.1    - paper Step 5

   These tools are NOT included in this repository.  They are governed
   by their own license, which restricts redistribution.  See
   [THIRD_PARTY.md](THIRD_PARTY.md) for how to obtain them, the citation
   requirements, and the data flow between their outputs and CShell's
   Python pipeline.

3. If you use CShell in a publication, please cite our paper [PLEASE
   FILL: paper citation once accepted] **and** the two Liu et al. papers
   listed in `THIRD_PARTY.md`.


--------------------------------------------------------------------------
II. INSTALLATION
--------------------------------------------------------------------------
**CShell is locked to the original PCshell notebook author's environment**
(Anaconda 2019.10; Python 3.7.4 + numpy 1.17.2 + pandas 0.25.1 +
networkx 2.3 + sympy 1.4).  This is **not** an arbitrary lower bound -
sympy 1.4's `Matrix(...).rref()` on a float64 matrix correctly zeroes
dependent rows, and newer sympy versions (>= ~1.5) produce numerically
tiny residues that break `cshell/step8_cp_eff.py` and yield NaN energy.
See `legacy_repro/REPORT.md` for the full investigation.

Recommended setup, via conda:

    git clone https://github.com/gaSuis/paper.git
    cd paper
    conda env create -f environment.yml
    conda activate cshell201910

(`requirements.txt` is also provided for pip users, but you must
ensure Python 3.7.4 is in use before installing.)

No compilation step is needed for the Python pipeline.  The third-party
Step 1 / Step 5 binaries from Yang-Yu Liu are only needed when you want
to **rebuild** the global classification or the per-disease
control-centrality files; for the asthma example shipped with this repo
those files are already pre-computed and stored under `data/`.


--------------------------------------------------------------------------
III. USAGE
--------------------------------------------------------------------------

### 1. Quick start - run the bundled asthma example

```bash
conda activate cshell201910
python scripts/run_pipeline.py \
    --network         data/network/Directed_network.txt \
    --disease-core    data/examples/asthma/disease_core.txt \
    --indispensable   data/network/Indispensable.txt \
    --extraordinary   data/network/Extraordinary.txt \
    --regulation-cc   data/examples/asthma/regulation_subgraph.txt.ControlCentrality \
    --observation-cc  data/examples/asthma/observation_subgraph.txt.ControlCentrality \
    --out-dir         results/asthma \
    --seed            0
```

Reference output for the asthma example is in
`data/examples/asthma/expected_output/`.  These reference files come from
the original notebooks' single unseeded run; they are kept for
inspection only.  Your output **will not be byte-identical** to those
files even with `--seed`, because the original notebooks did not use
deterministic RNG.  What you get with `--seed` is **bit-stable across
multiple runs of this pipeline on the same machine** (i.e. you can
re-run with the same seed and verify the result is reproducible).

A note on **CP-EFF**: Step 8 is a random search that records a gene set
only when the randomly-perturbed energy is *lower* than both the
previous random sample and the unperturbed original energy.  For the
asthma example the unperturbed original energy is already close to
optimal, so most random seeds produce **0 CP-EFF genes** - this is a
property of the published algorithm, not a bug.  Increase
`--n-random-step8` (default 100) to sample more candidates, or accept
the empty CP-EFF as a valid algorithmic output.


### 2. Inputs

| CLI flag           | File                                                  | Origin                            |
|--------------------|--------------------------------------------------------|-----------------------------------|
| `--network`        | 2-column directed edgelist                            | user-supplied                     |
| `--disease-core`   | 1 gene per line                                       | user-supplied                     |
| `--indispensable`  | 1 gene per line, `Indispensable.txt`                  | paper Step 2 from `.nodetype`     |
| `--extraordinary`  | 1 gene per line, `Extraordinary.txt`                  | paper Step 2 from `.nodetype`     |
| `--regulation-cc`  | `*.ControlCentrality` for the regulation subgraph     | third-party Step 5 (`CalControlCentrality`) |
| `--observation-cc` | `*.ControlCentrality` for the observation subgraph    | third-party Step 5                |

To rebuild `Indispensable.txt` / `Extraordinary.txt` for a new network,
or to rebuild the two `.ControlCentrality` files for a new disease core,
follow [THIRD_PARTY.md](THIRD_PARTY.md).


### 3. Outputs

Everything below is written under `--out-dir`:

```
cshell_genes.tsv         <-- unified result, two columns: gene  cp_class
                             cp_class is one of CP-IND / CP-EXT / CP-DRV /
                             CP-SEN / CP-EFF; one row per (gene, cp_class)

CP-IND.txt               per-class gene lists (one gene per line)
CP-EXT.txt
CP-DRV.txt               (= CP-DRV_part1 union CP-DRV_part2)
CP-DRV_part1.txt         intermediate from Step 6 (top 1% control subspace)
CP-DRV_part2.txt         intermediate from Step 7 (rank-fillers added on top of
                         CP-IND u CP-EXT)
CP-SEN.txt
CP-EFF.txt

regulation_subgraph.txt        <-- subgraph used by paper Step 5
observation_subgraph.txt       <-- subgraph used by paper Step 5
```

The `cshell_genes.tsv` file directly answers the
"which CP class selected gene X?" question without any post-processing:

```
gene     cp_class
ACSL1    CP-IND
...
HOXB8    CP-EFF
PAX4     CP-DRV
PAX4     CP-EFF
```

A gene that satisfies more than one CP definition appears on multiple
rows (one per CP class).


### 4. Run a single step

Use `scripts/run_step.py` for partial reruns (e.g. you changed only the
disease core but the network and its `.ControlCentrality` files are
unchanged).  See `python scripts/run_step.py -h` for the per-step
argument lists.


### 5. Pipeline overview

For a step-by-step description of the algorithm and how each step maps
to the paper, see [docs/pipeline_overview.md](docs/pipeline_overview.md).


--------------------------------------------------------------------------
IV. REPOSITORY LAYOUT
--------------------------------------------------------------------------
```
paper/
|-- README.md                        this file
|-- THIRD_PARTY.md                   third-party tools, licenses, citations
|-- LICENSE                          MIT (Python code only)
|-- environment.yml                  conda env spec (locked Anaconda 2019.10)
|-- requirements.txt                 pip-only fallback (versions match env yml)
|-- .gitignore
|-- cshell/                          our Python package
|   |-- __init__.py                  CP labels + legacy-name traceability tables
|   |-- io_utils.py
|   |-- step3_cp_ind_ext.py          CP-IND, CP-EXT
|   |-- step4_subgraphs.py           regulation / observation subgraph builder
|   |-- step6_cp_drv_sen.py          CP-DRV (part 1), CP-SEN
|   |-- step7_cp_drv_extra.py        CP-DRV (part 2)
|   |-- step8_cp_eff.py              CP-EFF
|   `-- pipeline.py                  orchestration + unified TSV writer
|-- scripts/
|   |-- run_pipeline.py              full pipeline CLI
|   `-- run_step.py                  single-step CLI
|-- legacy_repro/                    env-fidelity verification (reviewers)
|   |-- REPORT.md                    findings from rerunning the original code
|   |                                in the locked env
|   |-- inspect_notebooks.py         metadata + import inspector
|   |-- run_step8_legacy.py          notebook [5] verbatim port
|   |-- run_steps_3_6_7_legacy.py    notebooks [1] [3] [4] verbatim ports
|   |-- check_py37_compile.py        Python 3.7 syntax compatibility check
|   `-- legacy_*.txt                 raw outputs of the verbatim reruns
|-- data/
|   |-- network/
|   |   |-- Directed_network.txt
|   |   |-- Directed_network.txt.nodetype     (paper Step 1 output)
|   |   |-- Directed_network.txt.linktype     (paper Step 1 output)
|   |   |-- Indispensable.txt                 (paper Step 2 output)
|   |   |-- Extraordinary.txt                 (paper Step 2 output)
|   |   `-- Control_Signature_Result.txt      (paper Step 2 byproduct)
|   `-- examples/
|       `-- asthma/
|           |-- disease_core.txt
|           |-- regulation_subgraph.txt.ControlCentrality   (paper Step 5 output)
|           |-- observation_subgraph.txt.ControlCentrality  (paper Step 5 output)
|           `-- expected_output/                            (reference results)
`-- docs/
    `-- pipeline_overview.md
```

The two `*.ControlCentrality` files in `data/examples/asthma/` are the
output of the third-party `CalControlCentrality` tool when applied to
the `regulation_subgraph.txt` / `observation_subgraph.txt` produced by
paper Step 4 for the asthma core.  They are stored here for convenience
so the example can be reproduced **without** installing or running the
third-party C++ binaries.  See `THIRD_PARTY.md` for how to regenerate
them for a different disease core.
