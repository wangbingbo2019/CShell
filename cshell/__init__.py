"""CShell pipeline package.

Identify five classes of *control properties* (CPs) for a given disease core,
based on the network controllability framework of Liu et al. (Nature 2011) and
the control-centrality framework of Liu et al. (PLoS ONE 2012).

Five CP classes - paper terminology
-----------------------------------
The pipeline uses the upper paper labels consistently in all output files,
variable names, and the unified TSV's ``cp_class`` column.  Older
abbreviations from the Table S1 supplementary material (ISR_CS, ESR_CS,
CLS_CS, OLS_CS, MER_CS) and the very first version of the source code
(DCS_CS_1, DCS_CS_2, SOS_CS, IER_CS) are recorded below for traceability
only - they are NOT used as filenames in this repository.

  Paper label  |  Table S1 abbrev  |  Legacy code         |  Concept
  -------------|-------------------|----------------------|---------------------------
  CP-IND       |  ISR_CS           |  ISR_CS              |  Minimal driver set (indispensable)
  CP-EXT       |  ESR_CS           |  ESR_CS              |  Control signature (extraordinary)
  CP-DRV       |  CLS_CS           |  DCS_CS_1, DCS_CS_2  |  Control subspace
  CP-SEN       |  OLS_CS           |  SOS_CS              |  Control observability
  CP-EFF       |  MER_CS           |  IER_CS              |  Control energy
"""

__all__ = [
    "CP_CLASSES",
    "TABLE_S1_TO_CP",
    "LEGACY_TO_CP",
]

# canonical paper terminology, in pipeline order
CP_CLASSES = ("CP-IND", "CP-EXT", "CP-DRV", "CP-SEN", "CP-EFF")

# Table S1 abbreviation -> paper label (kept for traceability)
TABLE_S1_TO_CP = {
    "ISR_CS": "CP-IND",
    "ESR_CS": "CP-EXT",
    "CLS_CS": "CP-DRV",
    "OLS_CS": "CP-SEN",
    "MER_CS": "CP-EFF",
}

# legacy source-code variable -> paper label (kept for traceability)
LEGACY_TO_CP = {
    "ISR_CS":   "CP-IND",
    "ESR_CS":   "CP-EXT",
    "DCS_CS_1": "CP-DRV",
    "DCS_CS_2": "CP-DRV",
    "SOS_CS":   "CP-SEN",
    "IER_CS":   "CP-EFF",
}
