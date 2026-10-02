# Local reference QC source snapshots

See the [reference filtering guide](../../../docs/single-cell-reference-qc.md)
for inputs, runtime versions, commands, cohort-specific checks, and test
limits.

This directory contains the scDblFinder and SingleR source used in the
completed local Tabula Sapiens analyses, plus the final matrix-selection
code. The scientific scripts are unchanged. `source_manifest.json`
records their checksums. The source data and result matrices are external.

These are local Snakemake workflows, not Terra WDL workflows. There is no
Docker image build or CIBERSORTx calculation in this directory.
