# Full phenotype-PC outputs design

## Objective

Add full phenotype principal-component matrices to the prepare QTL pipeline without changing the existing Gavish–Donoho-selected PC files or downstream covariate behavior.

## Output contract

Each phenotype-PC calculation will continue to emit the existing selected-PC file:

- `<OutputPrefix>_phenotype_PCs<OutputSuffix>.tsv`

It will also emit a full rotated-PC matrix:

- `<OutputPrefix>_phenotype_PCs<OutputSuffix>.all.tsv`

For example, the INT and scaled outputs will be named `*_phenotype_PCs.INT.tsv`, `*_phenotype_PCs.INT.all.tsv`, `*_phenotype_PCs.scaled.tsv`, and `*_phenotype_PCs.scaled.all.tsv`.

The `.all.tsv` file will contain the same `ID` sample column and every available PC from the single PCA run. The existing file will retain the current Gavish–Donoho-selected columns and filename.

## Workflow and data flow

`scripts/common/calculate_PCs.R` will materialize both tables from the same PCA result. `workflows/common/calculate_phenotypePCs.wdl` will declare both files in the task output block and expose both through the workflow output block. This explicit declaration ensures Cromwell localizes/delocalizes the new full matrix instead of leaving it as an undeclared task artifact.

The expression, proteomics, splicing, and methylation workflows will expose the new full-PC outputs while continuing to pass the selected-PC outputs to existing covariate-merging calls. No existing output names or selected-PC consumers will change.

## Testing and documentation

Add regression coverage for the PC-generation behavior, checking that the full output preserves sample IDs and selected PC values and contains all available PCs. Add workflow-level checks that the `.all.tsv` files are declared and surfaced as outputs. Update the molecular-QTL and script documentation to describe both output classes and clarify that covariate merging continues to use the selected files.

## Error handling and compatibility

Both files must be written from the same successful PCA computation; a PCA failure should fail the task as it does today. The change is additive: existing callers that reference the current selected-PC filenames remain valid, while callers that need every phenotype PC can consume the new `.all.tsv` outputs.
