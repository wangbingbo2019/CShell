# CShell

CShell identifies control-property genes around a disease core in a directed
gene-regulatory network.  The pipeline reports five classes:

| Label | Meaning | Legacy name |
| --- | --- | --- |
| CP-IND | indispensable/minimal-driver related genes | ISR_CS |
| CP-EXT | extraordinary/control-signature genes | ESR_CS |
| CP-DRV | control-subspace genes | DCS_CS / CLS_CS |
| CP-SEN | control-observability genes | SOS_CS / OLS_CS |
| CP-EFF | control-energy reduction genes | IER_CS / MER_CS |

## Installation

Use the locked conda environment when possible:

```bash
conda env create -f environment.yml
conda activate cshell201910
```

The environment follows the original notebook runtime:

```text
python 3.7.4
numpy 1.17.2
pandas 0.25.1
networkx 2.3
sympy 1.4
```

`requirements.txt` is provided for pip-based installation, but Python 3.7
should be created first.

## Inputs

`scripts/run_pipeline.py` requires these files:

| Argument | Description |
| --- | --- |
| `--network` | two-column directed edgelist, one `source target` edge per line |
| `--disease-core` | disease core genes, one gene per line |
| `--indispensable` | global indispensable node list |
| `--extraordinary` | global extraordinary node list |
| `--regulation-cc` | `*.ControlCentrality` output for the regulation subgraph |
| `--observation-cc` | `*.ControlCentrality` output for the observation subgraph |
| `--out-dir` | directory for pipeline outputs |

The bundled network-level files, when distributed, are under `data/network/`.
Example disease-core and control-centrality files may be placed anywhere; pass
their paths explicitly on the command line.

## Example Command

```bash
python scripts/run_pipeline.py \
    --network         data/network/Directed_network.txt \
    --disease-core    data/asthma/disease_core.txt \
    --indispensable   data/network/Indispensable.txt \
    --extraordinary   data/network/Extraordinary.txt \
    --regulation-cc   data/asthma/regulation_subgraph.txt.ControlCentrality \
    --observation-cc  data/asthma/observation_subgraph.txt.ControlCentrality \
    --out-dir         results/asthma \
    --seed            0
```

The asthma sample data used by the paper consists of:

```text
disease_core.txt
regulation_subgraph.txt.ControlCentrality
observation_subgraph.txt.ControlCentrality
```

Optional reference outputs can be used for inspection, but they are snapshots
from one unseeded notebook run.  They are not guaranteed to be reproduced
byte-for-byte by a fresh run because Step 3 and Step 8 use random sampling.

## Outputs

The output directory contains:

```text
cshell_genes.tsv       gene-to-control-property table
CP-IND.txt             CP-IND genes
CP-EXT.txt             CP-EXT genes
CP-DRV.txt             union of CP-DRV_part1 and CP-DRV_part2
CP-DRV_part1.txt       top control-subspace genes
CP-DRV_part2.txt       rank-expansion control-subspace genes
CP-SEN.txt             control-observability genes
CP-EFF.txt             control-energy genes
regulation_subgraph.txt
observation_subgraph.txt
```

`cshell_genes.tsv` has two columns:

```text
gene    cp_class
```

A gene selected by more than one control property appears on multiple rows.

## Third-Party Tools

CShell uses outputs from two third-party C++ tools by Yang-Yu Liu.  They are
not redistributed here because they have separate license terms.

1. **ControllabilityAnalysis v1.0**
   Produces network-level controllability files such as `*.nodetype` and
   `*.linktype`.  These files are used to derive `Indispensable.txt` and
   `Extraordinary.txt`.

   Citation:
   Liu Y-Y, Slotine J-J, Barabasi A-L. Controllability of complex networks.
   Nature 473, 167-173 (2011).

2. **CalControlCentrality v1.1**
   Produces `*.ControlCentrality` files for the regulation and observation
   subgraphs.

   Citation:
   Liu Y-Y, Slotine J-J, Barabasi A-L. Control Centrality and Hierarchical
   Structure in Complex Networks. PLoS ONE 7(9): e44459 (2012).

See `THIRD_PARTY.md` for the expected data flow and citation details.

## Directory Layout

Core package layout:

```text
Cshell/
|-- README.md
|-- LICENSE
|-- environment.yml
|-- requirements.txt
|-- THIRD_PARTY.md
|-- cshell/
|   |-- io_utils.py
|   |-- pipeline.py
|   |-- step3_cp_ind_ext.py
|   |-- step4_subgraphs.py
|   |-- step6_cp_drv_sen.py
|   |-- step7_cp_drv_extra.py
|   `-- step8_cp_eff.py
|-- scripts/
|   |-- run_pipeline.py
|   `-- run_step.py
|-- data/
|   `-- network/
|       |-- Directed_network.txt
|       |-- Directed_network.txt.nodetype
|       |-- Directed_network.txt.linktype
|       |-- Indispensable.txt
|       |-- Extraordinary.txt
|       `-- Control_Signature_Result.txt
|-- docs/
|   `-- pipeline_overview.md
`-- supplementary/
```
