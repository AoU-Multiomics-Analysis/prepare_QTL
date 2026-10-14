# REGENIE exports

Add sample-row REGENIE files to the existing QTL preparation outputs. Cover
expression, protein, splicing, methylation, methylation cohort entry points,
cell type expression, and the phenotype-PC preparation entry point. Preserve
the existing BED and TensorQTL outputs.

A shared WDL 1.0 task reads a BED and optional covariates. It writes uncompressed
TSV phenotype batches (100 traits by default), a trait map with coordinates,
a batch table with the comma-separated phenotype selection, an optional
covariate TSV, a sample exclusion report, and a log. It uses FID=0 and the exact
sample ID as IID. Missing phenotypes become NA. Finite numeric text is retained.
Each trait gets a stable SHA-256-derived safe name; the map retains its source ID.
Duplicate source IDs, invalid values, and missing covariate cells cause an error.

INT and scaled branches use merged covariates when present and selected phenotype
PCs otherwise. Raw branches use additional covariates only. The sample set is the
intersection with the supplied covariates, in BED order; excluded IDs are reported.
The exporter does not normalize, impute, or residualize. Optional residualized BEDs
remain existing outputs; the shared exporter can be run on them separately.

The conversion keeps only one phenotype batch in memory and uses SQLite to check
duplicate IDs without retaining all feature IDs in memory. The task carries its
Python program in a literal heredoc, as in the existing phenotype-PC workflow.
Use the existing pinned REGENIE image, which supplies Python and REGENIE. No new
image is needed. Keep one source copy under tools and verify that it matches the
embedded program. Inputs remain File-typed until command rendering.

GitHub Actions validates WDL, exercises cloud-to-local command rendering, and
runs a synthetic two-step REGENIE scan using the exported files in the pinned
image. Do not build a local image or submit a Terra job. A local test or CI scan
does not establish that the complete workflow ran on Terra.
