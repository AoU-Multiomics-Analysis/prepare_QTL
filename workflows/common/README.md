# Common workflow helpers

These WDLs are shared across molecular modalities:

- `calculate_phenotypePCs.wdl` calculates phenotype PCs from a molecular BED and emits both Gavish-Donoho-selected and full rotated-PC TSV outputs.
- `MergeCovariates.wdl` combines the selected phenotype PCs with sample-level covariates in TensorQTL orientation.
- `ResidualizePhenotypes.wdl` residualizes and scales normalized phenotype BEDs.

Modality workflows import these helpers with paths relative to their own directory, for example `../common/calculate_phenotypePCs.wdl`.
