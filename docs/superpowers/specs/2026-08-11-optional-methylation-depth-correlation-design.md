# Optional methylation-depth correlation design

## Goal

Allow the PacBio methylation preprocessing workflows to skip the per-CpG Spearman correlation between methylation abundance and normalized sequencing depth while preserving the existing output schema and all other QC behavior.

## Scope

The option applies only to `coverage_methylation_spearman_rho` and `n_samples_coverage_methylation_correlation`. It does not change minimum-coverage filtering, required extreme-coverage filtering, local CpG-correlation clustering, sample-connectivity filtering, site filters, or phenotype generation.

## Design

Add `ComputeCoverageMethylationCorrelation` as a Boolean WDL input with default `true` to the current cohort workflow and the legacy manifest/array cohort entry points. Propagate the value to the Rust chromosome merger and the R chromosome merger. The Rust CLI will expose a `--skip-coverage-methylation-correlation` switch; the R CLI will expose the equivalent `--SkipCoverageMethylationCorrelation` switch.

When computation is enabled, outputs are unchanged. When disabled, the cohort metadata retains the same two columns with `coverage_methylation_spearman_rho` set to `NA` and `n_samples_coverage_methylation_correlation` set to `0`. The implementation will avoid collecting correlation vectors and avoid calculating normalized coverage solely for this diagnostic when disabled.

## Testing

Add Rust unit coverage for the enabled and disabled correlation-summary behavior, plus command-line parsing coverage for the skip flag. Run Rust tests and validate every modified WDL parses successfully. Inspect the generated source diff to confirm no unrelated QC behavior changed.
